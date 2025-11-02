#!/bin/bash
set -e
set -u
set -o pipefail

# Filter samples & calculate PC

DIR="/data/projects/endometriosis/tensorQTL_final/final_20250210/results/imp_130_final"

if [ ! -d "$DIR" ]; then
  mkdir -p "$DIR"
fi

# Keep only children samples

plink2 \
  --bfile /data/projects/endometriosis/qc_final/imputed_159_samples_rsq09_maf001_hwe005_chr1-22_sorted_wo_multiallelic \
  --keep /data/projects/endometriosis/tensorQTL_final/final_20250210/results/plink_130_ids.tsv \
  --make-bed \
  --out ${DIR}/imp_130_final


# Calculate PCs

plink \
  --bfile ${DIR}/imp_130_final \
  --indep-pairwise 50 5 0.2 \
  --out ${DIR}/imp_130_final

plink2 \
  --bfile ${DIR}/imp_130_final \
  --extract ${DIR}/imp_130_final.prune.in \
  --pca \
  --out ${DIR}/imp_130_final