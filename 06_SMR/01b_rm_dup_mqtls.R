
setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_nom")

library(tidyverse)

d <- read_tsv("tensorqtl_130_5e-08_20250210.tsv")

uniq_snps <- unique(d$SNP)
chr <- str_split_i(uniq_snps, ":", 1)
pos <- str_split_i(uniq_snps, ":", 2)

chrpos <- paste0(chr, ":", pos)

rm_dups <- uniq_snps[duplicated(chrpos)]

# d <- d %>%
#     filter(str_length(A1) <= 1) %>%
#     filter(str_length(A2) <= 1) %>%
#     mutate(qtl = paste0(SNP, "_", Probe)) %>%
#     filter(!duplicated(qtl)) %>%
#     select(-qtl) %>% 
#     mutate(SNP = paste0(Chr, ":", BP)) %>% 
#     mutate(qtl = paste0(SNP, "_", Probe)) %>% 
#     filter(!duplicated(qtl)) %>% 
#     select(-qtl)


d <- d %>% 
    filter(!SNP %in% rm_dups) %>% 
    mutate(SNP = paste0(Chr, ":", BP)) %>% 
    mutate(qtl = paste0(SNP, "_", Probe)) %>% 
    filter(!duplicated(qtl)) %>%
    select(-qtl)

#write_tsv(d, "./endometriosis_20240521/endo_pval5e-08_smrquery.rm_dup_mqtls.chrpos.tsv")
write_tsv(d, "tensorqtl_130_5e-08_20250210_rm_dup_mqtls_chrpos.tsv")
