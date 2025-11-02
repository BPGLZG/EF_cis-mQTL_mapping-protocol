
setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210")

library(tidyverse)

# Use already prepared cov file to filter the 130 samples from the 159 file
cov <- read_tsv("data/covs_130_tensor.txt")
cov_samples <- colnames(cov)[-1]

# Missing FID in the names. Using original plink .fam file
fam <- read_table("/data/projects/endometriosis/qc_final/imputed_159_samples_rsq09_maf001_hwe005_chr1-22_sorted_wo_multiallelic.fam",
                col_names = FALSE)

# Filter 130 rows and select FID and IID (for the file need in the --keep flag in PLINK)

ids <- fam %>%
  select(X1, X2) %>%
  filter(X2 %in% cov_samples)

# Save results
write_tsv(ids, "results/plink_130_ids.tsv", col_names = FALSE)
