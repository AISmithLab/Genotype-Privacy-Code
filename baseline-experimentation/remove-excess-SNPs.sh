#!/bin/bash
max_layer_removed=$1
# ==============================================================
# Adjustable parameters for experimenting with different samples
sample="HG00108"
source_snp="rs9315973"
snp_chromosome="chr13"
snp_query="13:43501356"
num_experiments=3
# ==============================================================
to_remove="./not_directly_genotyped_${source_snp}_neighbors.txt"
DATA_DIR="../experiment_data"
mkdir -p "$DATA_DIR"
RESULTS_DIR="./baseline_01_${source_snp}"
mkdir -p "$RESULTS_DIR"
RESULTS_PATH="${RESULTS_DIR}/dn_${sample}_imputation_results_remove_degree_${max_layer_removed}.txt"
for exp_num in $(seq 1 "$num_experiments"); do
    (
        bcftools view -e \
        "ID=@${to_remove}" \
        -O z \
        -o "${DATA_DIR}/${sample}-${snp_chromosome}-exp${exp_num}-${source_snp}-to-impute.vcf.gz" \
        "../chromosomes/${sample}-${snp_chromosome}.vcf.gz"
    ) &
done

wait

for exp_num in $(seq 1 "$num_experiments"); do
    (
        java -jar beagle.29Oct24.c8e.jar \
        gt="${DATA_DIR}/${sample}-${snp_chromosome}-exp${exp_num}-${source_snp}-to-impute.vcf.gz" \
        gp=true \
        ref="../chromosomes/no-${sample}-${snp_chromosome}.vcf.gz" \
        out="${DATA_DIR}/${sample}-${snp_chromosome}-exp${exp_num}-imputed" \
        && tabix -p vcf "${DATA_DIR}/${sample}-${snp_chromosome}-exp${exp_num}-imputed.vcf.gz" \
        && bcftools view -r "${snp_query}" \
        "${DATA_DIR}/${sample}-${snp_chromosome}-exp${exp_num}-imputed.vcf.gz" | grep -v '^##' >> "$RESULTS_PATH"
    ) &
done

wait
