# SMR_final - 08_filter_smr.R (.msmr to .tsv)
# February 11, 2025

# a) Breast cancer (ieu-a-1126)
# b) Breast cancer (GCST004988)
# c) Breast cancer (GCST010098)
# d) Cervical cancer (GCST004833)
# e) Endometrial cancer (GCST006464)
# f) Endometriosis (GCST90269970)
# g) Endometriosis (GCST90205183)
# h) Epithelial ovarian cancer (GCST004462)***
# i) Epithelial ovarian cancer (GCST90244167)
# j) Ovarian cancer (ieu-a-1120)
# k) Uterine fibroids (GCST009158)

# Set working directory
setwd("/data/projects/endometriosis/smr_final/smr_130_5_9/")

# Load libraries
library(data.table)
library(tidyverse)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#!/usr/bin/env Rscript

# Read SMR results
#-------------------------------------------------------------------------------
# a) Breast cancer (ieu-a-1126)
ressmr_a <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/ieu-a-1126/ieu-a-1126.msmr")

ressmr_a <- ressmr_a %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_a <- ressmr_a %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_a, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/ieu-a-1126/smr_ieu-a-1126_filtered.tsv")

# Load filtered smr results
sigres_a <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/ieu-a-1126/smr_ieu-a-1126_filtered.tsv")


# Write CpG info
anno_a <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_a <- anno_a %>% 
  filter(Name %in% sigres_a$probeID)

write_tsv(anno_sigres_a, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/ieu-a-1126/smr_ieu-a-1126_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# b) Breast cancer (GCST004988)
#ressmr_b <- read_tsv("/data/projects/endometriosis/smr_final/20250120_smr_145_5_20/results/breast_cancer/GCST004988/GCST004988.msmr")
ressmr_b <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST004988/GCST004988.msmr")

ressmr_b <- ressmr_b %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_b <- ressmr_b %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_b, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST004988/smr_GCST004988_filtered.tsv")

# Load filtered smr results
sigres_b <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST004988/smr_GCST004988_filtered.tsv")


# Write CpG info
anno_b <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_b <- anno_b %>% 
  filter(Name %in% sigres_b$probeID)

write_tsv(anno_sigres_b, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST004988/smr_GCST004988_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# c) Breast cancer (GCST010098)
ressmr_c <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/GCST010098.msmr")

ressmr_c <- ressmr_c %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_c <- ressmr_c %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_c, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/smr_GCST010098_filtered.tsv")

# Load filtered smr results
sigres_c <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/smr_GCST010098_filtered.tsv")


# Write CpG info
anno_c <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_c <- anno_c %>% 
  filter(Name %in% sigres_c$probeID)

write_tsv(anno_sigres_c, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/smr_GCST010098_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# d) Cervical cancer (GCST004833)
ressmr_d <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/cervical_cancer/GCST004833/GCST004833.msmr")

ressmr_d <- ressmr_d %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_d <- ressmr_d %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_d, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/cervical_cancer/GCST004833/smr_GCST004833_filtered.tsv")

# Load filtered smr results
sigres_d <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/cervical_cancer/GCST004833/smr_GCST004833_filtered.tsv")


# Write CpG info
anno_d <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_d <- anno_d %>% 
  filter(Name %in% sigres_d$probeID)

write_tsv(anno_sigres_d, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/cervical_cancer/GCST004833/smr_GCST004833_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# e) Endometrial cancer (GCST006464)
ressmr_e <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/GCST006464.msmr")

ressmr_e <- ressmr_e %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_e <- ressmr_e %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_e, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/smr_GCST006464_filtered.tsv")

# Load filtered smr results
sigres_e <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/smr_GCST006464_filtered.tsv")


# Write CpG info
anno_e <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_e <- anno_e %>% 
  filter(Name %in% sigres_e$probeID)

write_tsv(anno_sigres_e, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/smr_GCST006464_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# f) Endometriosis (GCST90269970)
ressmr_f <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/GCST90269970.msmr")

ressmr_f <- ressmr_f %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_f <- ressmr_f %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_f, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/smr_GCST90269970_filtered.tsv")

# Load filtered smr results
sigres_f <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/smr_GCST90269970_filtered.tsv")


# Write CpG info
anno_f <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_f <- anno_f %>% 
  filter(Name %in% sigres_f$probeID)

write_tsv(anno_sigres_f, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/smr_GCST90269970_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# g) Endometriosis (GCST90205183)
ressmr_g <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/GCST90205183.msmr")

ressmr_g <- ressmr_g %>% 
    mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_g <- ressmr_g %>% 
    filter(p_SMR_multi_bonf < 0.05) %>% 
    filter(p_HEIDI > 0.05) %>% 
    arrange(p_SMR_multi)

write_tsv(sigres_g, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/smr_GCST90205183_filtered.tsv")

# Load filtered smr results
sigres_g <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/smr_GCST90205183_filtered.tsv")


# Write CpG info
anno_g <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
    as_tibble()
anno_sigres_g <- anno_g %>% 
    filter(Name %in% sigres_g$probeID)

write_tsv(anno_sigres_g, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/smr_GCST90205183_cpgs_info.tsv")

# ------------------------------------------------------------------------------
# h) Epithelial ovarian cancer (GCST004462)
ressmr_h <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST004462/GCST004462.msmr")

ressmr_h <- ressmr_h %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_h <- ressmr_h %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_h, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST004462/smr_GCST004462_filtered.tsv")

# Load filtered smr results
sigres_h <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST004462/smr_GCST004462_filtered.tsv")


# Write CpG info
anno_h <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_h <- anno_h %>% 
  filter(Name %in% sigres_h$probeID)

write_tsv(anno_sigres_h, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST004462/smr_GCST004462_cpgs_info.tsv")

# ------------------------------------------------------------------------------
# i) Epithelial ovarian cancer (GCST90244167)
ressmr_i <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST90244167/GCST90244167.msmr")

ressmr_i <- ressmr_i %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_i <- ressmr_i %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_i, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST90244167/smr_GCST90244167_filtered.tsv")

# Load filtered smr results
sigres_i <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST90244167/smr_GCST90244167_filtered.tsv")


# Write CpG info
anno_i <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_i <- anno_i %>% 
  filter(Name %in% sigres_i$probeID)

write_tsv(anno_sigres_i, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/epithelial_oc/GCST90244167/smr_GCST90244167_cpgs_info.tsv")

# ------------------------------------------------------------------------------
# j) Ovarian cancer (ieu-a-1120)
ressmr_j <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/ieu-a-1120.msmr")

ressmr_j <- ressmr_j %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_j <- ressmr_j %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_j, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/smr_ieu-a-1120_fitered.tsv")

# Load filtered smr results
sigres_j <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/smr_ieu-a-1120_fitered.tsv")


# Write CpG info
anno_j <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_j <- anno_j %>% 
  filter(Name %in% sigres_j$probeID)

write_tsv(anno_sigres_j, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/smr_ieu-a-1120_cpgs_info.tsv")

#-------------------------------------------------------------------------------
# k) Uterine fibroids (GCST009158)
ressmr_k <- read_tsv("/data/projects/endometriosis/smr_final/smr_130_5_9/results/uterine_fibroids/GCST009158/GCST009158.msmr")

ressmr_k <- ressmr_k %>% 
  mutate(p_SMR_multi_bonf = p.adjust(p_SMR_multi, method = "bonferroni"))

sigres_k <- ressmr_k %>% 
  filter(p_SMR_multi_bonf < 0.05) %>% 
  filter(p_HEIDI > 0.05) %>% 
  arrange(p_SMR_multi)

write_tsv(sigres_k, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/uterine_fibroids/GCST009158/smr_GCST009158_filtered.tsv")

# Load filtered smr results
sigres_k <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/uterine_fibroids/GCST009158/smr_GCST009158_filtered.tsv")


# Write CpG info
anno_k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble()
anno_sigres_k <- anno_k %>% 
  filter(Name %in% sigres_k$probeID)

write_tsv(anno_sigres_k, "/data/projects/endometriosis/smr_final/smr_130_5_9/results/uterine_fibroids/GCST009158/smr_GCST009158_cpgs_info.tsv")

