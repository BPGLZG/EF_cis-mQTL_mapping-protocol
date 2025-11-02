### cis-mQTL calculation - EFA
# February 10, 2025

# Set working directory
setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210")

# Load required libraries
library(tidyverse)
library(dplyr)

# Load methylation values
met <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_rnt_qced.bed")

dim(met)    #626262 CpGs + 149 samples

# Load covariates
cov <- read_tsv("data/covs_130_tensor.txt",
                col_types = cols(.default = col_character()))

dim(cov)    #4 covariates + 131 samples

# Load methylation PCs
mpcs <- read_tsv("data/residualized_11_mPC.transposed.tsv",
                 col_types = cols(.default = col_character()))

dim(mpcs)   #9 mPCs + 131 samples

# Load sample names to keep and reorder only the prepared samples
final_samples <- colnames(cov)[-1]

# Select only the prepared samples
met <- met %>%
  select(all_of(c(colnames(met)[1:4], final_samples)))

dim(met)   #626262 CpGs + 134 sampes

# Load new genotype PCs
pcs <- read_tsv("results/imp_130_final/imp_130_final.eigenvec", trim_ws = TRUE)
pcs <- t(as.matrix(pcs[, -1]))
colnames(pcs) <- pcs[1, ]
pcs <- pcs[-1, ]
pcs <- as_tibble(pcs, rownames = "id")
pcs <- pcs[1:5, ]

# Bind covariates, PCs and mPCs
colnames(cov)[1] <- "id"
finalcov <- bind_rows(cov, pcs, mpcs)

dim(finalcov)   # 18 covs (4covs + 5gPCs + 9mPCs) + 131 samples

write_tsv(met, "results/tensorqtl_130_met.bed.gz") #130 samples + #Chr, start, end, ID
write_tsv(finalcov, "results/tensorqtl_covs.tsv")  #4 covs + 5gPCs + 9mPCs
