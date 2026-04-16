#!/bin/bash

(
    ./plink_linux_x86_64/plink \
    --vcf ../chromosomes/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --make-bed --out ../plink_outputs/1kGP-chr11
)

wait

for ld_threshold in 0.05; do
    (
        ./plink_linux_x86_64/plink \
            --bfile ../plink_outputs/1kGP-chr11 \
            --r2 --ld-window 10000 \
            --ld-window-kb 10000 \
            --ld-window-r2 "${ld_threshold}" \
            --with-freqs --out "../plink_outputs/1kGP-chr11-LD-min-${ld_threshold}"
        mv "../plink_outputs/1kGP-chr11-LD-min-${ld_threshold}.ld" "../plink_outputs/1kGP-chr11-LD-min-${ld_threshold}.txt"
    )
done
