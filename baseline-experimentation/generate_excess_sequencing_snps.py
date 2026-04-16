from cyvcf2 import VCF
import pandas as pd
import argparse
from tqdm import tqdm

######################################################
# DATA_FILEPATHS contains adjustable parameters
######################################################
DATA_FILEPATHS = {
    'manifest': '../illumina_manifests/GSA-manifest-build37.csv',
    'target_vcf': '../chromosomes/HG00096-chr13.vcf.gz',
    'source_snp': 'rs9315973',
    'source_snp_id': 'rs9315973',
    'source_snp_chromosome': '13',
    'min_ld_association_strength': '0.05'
}

def init_argparser():
    parser = argparse.ArgumentParser(
        description="Arguments parser for removing SNPs not found in Illumina manifest."
    )
    parser.add_argument(
        "-n", "--neighbor-max-degree",
        type=int,
        default=10,
        help="Remove at maximum n-th degree neighboring SNPs (default: %(default)d)"
    )
    parser.add_argument(
        "-a", "--attack",
        action="store_true",
        help="Enable imputation attack by removing neighboring SNPs to target SNP."
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    command_line_args = init_argparser()
    
    df = pd.read_csv(DATA_FILEPATHS['manifest'])
    df = df[df['Chr'].notna()]
    df = df[['IlmnID', 'Name', 'SNP', 'Chr', 'MapInfo']]
    df['MapInfo'] = df['MapInfo'].astype(int)
    df['Chr'] = df['Chr'].astype(str)
    
    vcf_filepath = DATA_FILEPATHS['target_vcf']
    vcf = VCF(vcf_filepath, gts012=True)

    meta_info = {'CHROM':[], 'POS':[], 'REF':[], 'ALT':[], 'ID': [], 'AlleleFreq': []}
    for v in tqdm(vcf):
        meta_info['CHROM'].append(v.CHROM)
        meta_info['POS'].append(v.POS)
        meta_info['REF'].append(v.REF)
        meta_info['ALT'].append(v.ALT[0])
        meta_info['ID'].append(v.ID)
        meta_info['AlleleFreq'].append(v.INFO['AF'])
        
    meta = pd.DataFrame(meta_info)
    neighbor_max_degree = command_line_args.neighbor_max_degree
    target_neighbors = []
    with open(f"./{DATA_FILEPATHS['source_snp']}-neighbors-by-degree/{neighbor_max_degree}-degree-neighbors-r2-{DATA_FILEPATHS['min_ld_association_strength']}.txt", 'r') as f:
        target_neighbors.extend([line.strip() for line in f.readlines()])
        
    genotyped_chr17 = df[df['Chr'] == DATA_FILEPATHS['source_snp_chromosome']]
    print(f"{genotyped_chr17.shape[0]} SNPs are genotyped by Illumina according to manifest")
    not_directly_genotyped = meta[~meta['POS'].isin(genotyped_chr17['MapInfo'])]
    directly_genotyped = meta[meta['POS'].isin(genotyped_chr17['MapInfo'])]
    target_neighbors_typed = directly_genotyped[directly_genotyped['ID'].isin(pd.Series(target_neighbors))]
    print(f"{target_neighbors_typed.shape[0]} SNPs in {DATA_FILEPATHS['source_snp']} {neighbor_max_degree}th degree neighbors are genotyped in 1kGP reference")
    print(f"{not_directly_genotyped.shape[0]} from 1kGP reference are not genotyped by Illumina")
    print(f"{directly_genotyped.shape[0]} SNPs from manifest are found in 1kGP reference")
    
    with open(f"not_directly_genotyped_{DATA_FILEPATHS['source_snp']}_neighbors.txt", 'w') as f:
        if command_line_args.attack:
            not_typed_to_remove = pd.concat([not_directly_genotyped, target_neighbors_typed], ignore_index=True)
        else:
            not_typed_to_remove = not_directly_genotyped
        counter = 0
        for snp in not_typed_to_remove['ID']:
            f.write(f'{snp}\n')
            counter += 1
        if DATA_FILEPATHS['source_snp_id'] not in set(not_typed_to_remove['ID'].to_list()):
            f.write(f"{DATA_FILEPATHS['source_snp_id']}\n")
            counter += 1
        print('='*30)
        print(f'{counter} SNPs not directly genotyped according to manifest and removed')
        print('='*30)
        
    with open(f"directly_genotyped_{DATA_FILEPATHS['source_snp']}_neighbors.txt", 'w') as f:
        for snp in directly_genotyped['ID']:
            if snp != DATA_FILEPATHS['source_snp_id']:
                f.write(f"{snp}\n")
