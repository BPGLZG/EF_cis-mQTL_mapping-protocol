# SMR_final - 07_run_smr.sh (output = .msmr file)
# February 03, 2025

# a) Breast cancer (ieu-a-1126)
# b) Breast cancer (GCST004988)
# c) Breast cancer (GCST010098)
# d) Cervical cancer (GCST004833)
# e) Endometrial cancer (GCST006464)
# f) Endometriosis (GCST90269970)
# g) Endometriosis (GCST90205183)
# h) Epithelial ovarian cancer (GCST004462)
# i) Epithelial ovarian cancer (GCST90244167)
# j) Ovarian cancer (ieu-a-1120)
# k) uterine fibroids (GCST009158)

# Working directory
"/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/"
/data/projects/endometriosis/smr_final_

#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Convert SMR Query to BESD (to do only once)
# Generates .esi, .epi, and .besd files

smr \
    --qfile "/data/projects/endometriosis/smr_final/smr_130_5_9/data/endo_pval5e-08_smrquery.tsv" \
    --make-besd \
    --out "/data/projects/endometriosis/smr_final/smr_130_5_9/data/endo_pval5e-08"

#-------------------------------------------------------------------------------
# Run SMR

REF="/data/projects/inma_hrc_imputation/imputed/chrAll"

BEQTL="/data/projects/endometriosis/smr_final/smr_130_5_9/data/endo_pval5e-08"

# GWAS datasets
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/ieu-a-1126_breast_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST004988_breast_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST010098_breast_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST004833_cervical_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST006464_new_endometrial_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST90269970_endometriosis.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST90205183_endometriosis.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST004462_epithelial_oc.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST90244167_epithelial_oc.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/ieu-a-1120_new_ovarian_cancer.ma"
GWAS="/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/data/gwas/GCST009158_uterine_fibroids.ma"

# Output files
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/ieu-a-1126/ieu-a-1126"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST004988/GCST004988"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/GCST010098"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/cervical_cancer/GCST004833/GCST004833"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/GCST006464"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/GCST90269970"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/GCST90205183"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST004462/GCST004462"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST90244167/GCST90244167"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/ieu-a-1120"
OUTNAME="/data/projects/endometriosis/smr_final/smr_130_5_9/results/uterine_fibroids/GCST009158/GCST009158"

#SMR command
smr --bfile $REF \
    --gwas-summary $GWAS \
    --beqtl-summary $BEQTL \
    --smr-multi \
    --thread-num 16 \
    --out ${OUTNAME} \
    > ${OUTNAME}.stdout.txt \
    2> ${OUTNAME}.stderr.txt
    