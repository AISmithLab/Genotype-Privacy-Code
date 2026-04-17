import matplotlib.pyplot as plt
import numpy as np
import glob

if __name__ == '__main__':
    
    ##########################################################
    # Adjustable parameter for plotting different SNPs and 
    # different baseline depth
    SNP_ID = 'rs78540526'
    MAX_DEPTH = 15
    ##########################################################
    RES_DIR = f'../baseline-experimentation/baseline_01_{SNP_ID}'
    DEPTH_RANGE = list(range(1,16))
    res_files_1 = sorted(glob.glob(f'{RES_DIR}/dn_*.txt'), key=lambda x: int(x.split('_')[-1].split('.')[0]))
    res_files_2 = sorted(glob.glob(f'{RES_DIR}/nd_*.txt'), key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharex=True)
    labels = ['Risk Dosage (mean +/- 1 std)', 'Normal Dosage (mean +/- 1 std)']
    colors = ['#d81159', '#0496ff']
    
    for fileset in [res_files_1, res_files_2]:
        counter = 0
        increment = len(fileset) // MAX_DEPTH
        num_variants_mean = []
        num_variants_std = []
        genotype_calling_accuracy = []
        while counter < len(fileset):
            expected_num_variants_res = []
            called_genotypes = []
            working_file = fileset[-1]
            for res_f in fileset[counter:counter+increment]:
                with open(res_f, 'r') as f:
                    # each individual file main contain duplicate results, and
                    # imputation of the same state is deterministic, so only keep 1
                    genotype_stats = f.readlines()[1].split()[-1].split(':')
                    called_genotypes.append(genotype_stats[0])
                    expected_num_variants = genotype_stats[1]
                    expected_num_variants_res.append(float(expected_num_variants))
            num_variants_mean.append(np.mean(expected_num_variants_res))
            num_variants_std.append(np.std(expected_num_variants_res))
            if 'dn' in working_file:
                genotype_calling_accuracy.append(len([i for i in called_genotypes if (i == '1|1' or i == '1|0' or i == '0|1')]) / len(called_genotypes))
            else:
                genotype_calling_accuracy.append(len([i for i in called_genotypes if (i == '0|0')]) / len(called_genotypes))
            counter += increment
        
        # Clip error bars to y-axis limits (0, 2)
        num_variants_mean = np.array(num_variants_mean)
        num_variants_std = np.array(num_variants_std)
        lower_err = np.minimum(num_variants_std, num_variants_mean - 0)  # clip at y=0
        upper_err = np.minimum(num_variants_std, 1.999 - num_variants_mean)  # clip at y=2
        
        if 'dn' in fileset[0]:
            ax1.errorbar(DEPTH_RANGE, num_variants_mean, [lower_err, upper_err], fmt='-o', capsize=5, capthick=1.5, linewidth=2, color=colors[0], label=labels[0])
            ax2.plot(DEPTH_RANGE, genotype_calling_accuracy, marker='o', linewidth=2, color=colors[0], label='Imputation Accuracy for Risk Samples')
        else:
            ax1.errorbar(DEPTH_RANGE, num_variants_mean, [lower_err, upper_err], fmt='-o', capsize=5, capthick=1.5, linewidth=2, color=colors[1], label=labels[1])
            ax2.plot(DEPTH_RANGE, genotype_calling_accuracy, marker='o', linewidth=2, color=colors[1], label='Imputation Accuracy for Control Samples')

    ax1.set_xticks(list(range(1,16)))
    ax2.set_xticks(list(range(1,16)))
    ax1.tick_params(axis='both', labelsize=18)
    ax2.tick_params(axis='both', labelsize=18)
    ax2.set_xlabel('Depth of Neighbors Removed', fontsize=20)
    ax1.set_ylim(-0.01, 2.01)
    ax2.set_ylim(-0.01, 1.01)
    ax1.set_ylabel('Allele Dosage', fontsize=20)
    ax2.set_ylabel('Imputation Accuracy', fontsize=20)
    ax1.legend(loc='best', fontsize=20)
    ax2.legend(loc='best', fontsize=20)
    ax1.grid(True)
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig(f'./output_figures/baseline-remove-onebit-attack-{SNP_ID}.png', dpi=300, bbox_inches='tight')
