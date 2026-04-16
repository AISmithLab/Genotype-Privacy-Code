#!/bin/bash
for threshold in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
    (
        python3 generate_excess_sequencing_snps.py -a -n "${threshold}" \
        && ./remove-excess-SNPs.sh "${threshold}"
    )
done