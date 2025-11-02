setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210/")

library(tidyverse)

# Load methylation Rank Based INT values
met <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_rnt_qced.bed",
                num_threads = 16)

met_t <- met %>%
  select(-`#Chr`, -start, -end) %>%
  t() %>%
  as.data.frame()

colnames(met_t) <- as.character(as.vector(met_t[1, ]))
met_t <- met_t[-1, ]
sample_id <- row.names(met_t)
met_t <- apply(met_t, 2, as.numeric)
rownames(met_t) <- sample_id
rm(sample_id)

met_t <- as_tibble(met_t, rownames = "IID")

# Load phenotype data
# pheno <- read_tsv("/data/projects/Placenta/EPIC/PACEanalysis/finalsamples_FID_IID_Basename.txt")

# Load covariates
pheno <- read_tsv("data/covs_130_tensor.txt")
pheno <- t(pheno)
colnames(pheno) <- pheno[1, ]
pheno <- pheno[-1, ]
pheno <- as_tibble(as.data.frame(pheno),
                   rownames = "IID")

# Ensure class
pheno$Age <- as.numeric(pheno$Age)
pheno$Batch <- as.factor(as.integer(pheno$Batch))
pheno$Group <- as.factor(as.integer(pheno$Group))
pheno$Phase <- as.factor(as.integer(pheno$Phase))

# Load genotype PCs
pcs <- read_tsv("results/imp_130_final/imp_130_final.eigenvec") %>%
  select(IID, paste0("PC", 1:5))


# Merge data
m <- inner_join(
  x = pheno,
  y = pcs,
  by = join_by("IID")
)

m <- inner_join(
  x = m,
  y = met_t,
  by = join_by("IID")
)


# Calculate residuals
m_met <- m %>%
  select(-colnames(pheno), -colnames(pcs)) %>%
  as.matrix()

fit <- lm(m_met ~ Age + Batch + Group + Phase + PC1 + PC2 + PC3 + PC4 + PC5, data = m)
residuals <- as.matrix(fit$residuals)

# Perform PCA on the residuals (mPCs)
pca <- prcomp(residuals)
colnames(pca$x) <- paste0("m", colnames(pca$x))

# Get variance explained per mPC
prop <- pca$sdev^2 / sum(pca$sdev^2) * 100
cumprop <- cumsum(pca$sdev^2) / sum(pca$sdev^2) * 100
pcname <- paste0("mPC", seq_along(prop))
pc <- data.frame(pcname, prop, cumprop)

# Select PCs explaining 80% of the accumulative variance,
# maximum of 20 (57 obs)
pc <- pc[cumprop <= 80, ]

pc <- pca$x[, colnames(pca$x) %in% pc$pcname]
max_pc <- if (ncol(pc) <= 20) {
  ncol(pc)
} else {
  20
}

pca20 <- pc[, 1:max_pc]
pca20 <- as.data.frame(pca20)
pca20 <- cbind("id" = m$IID, pca20)
write_tsv(pca20, "data/residualized_mPCs.tsv")

# Prepare transposed
pca20_t <- t(pca20) %>%
  as.data.frame()

write.table(pca20_t, "data/residualized_mPC.transposed.tsv",
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
