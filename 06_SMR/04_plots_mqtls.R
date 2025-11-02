#!/usr/bin/env Rscript

setwd("/data/projects/endometriosis/smr_final/smr_130_5_9")

#install.packages("rcartocolor")
library(rcartocolor)
library(tidyverse)
library(data.table)

# SMR results show that SNPs 2:11723110:A_G & 2:11731846:C_T are significant
# In PLINK files they are encoded as 2:11723110:G:A and 2:11731846:T:C

# SMR results show that SNPs 2:11723858:A_G, 11:30382899:C_T, 6:151767907:A_G, and 12:54387947:A_G are significant

# Genotype data
# system("plink --bfile ../data/plink/whole_maf05_hwe005_inmujeres --snp 2:11723110:G:A --recode --out ../data/plink/2_11723110_A_G")
# system("plink --bfile ../data/plink/whole_maf05_hwe005_inmujeres --snp 2:11731846:T:C --recode --out ../data/plink/2_11731846_C_T")

### Get the genotype data

# Genotype data 2:11723858:A_G    cg07314298
# system("plink --bfile /data/projects/endometriosis/smr_final/smr_130_5_9/data/sub_plink/imp_130_final --snp 2:11723858:A_G --recode --out /data/projects/endometriosis/smr_final/smr_130_5_9/data/plink/2:11723858:A_G)

# Genotype data 11:30382899:C_T   cg21393004
# system("plink --bfile /data/projects/endometriosis/smr_final/smr_130_5_9/data/sub_plink/imp_130_final --snp 11:30382899:C_T --recode --out /data/projects/endometriosis/smr_final/smr_130_5_9/data/plink/11:30382899:C_T)

# Genotype data 6:151767907:A_G   cg23088022
# system("plink --bfile /data/projects/endometriosis/smr_final/smr_130_5_9/data/sub_plink/imp_130_final --snp 6:151767907:A_G --recode --out /data/projects/endometriosis/smr_final/smr_130_5_9/data/plink/6:151767907:A_G)

# Genotype data 12:54387947:A_G   cg09605287
# system("plink --bfile /data/projects/endometriosis/smr_final/smr_130_5_9/data/sub_plink/imp_130_final --snp 12:54387947:A_G --recode --out /data/projects/endometriosis/smr_final/smr_130_5_9/data/plink/12:54387947:A_G)


### Format each SNP

# snp_1 - 2:11723858:A_G    cg07314298
snp1 <- read_delim("../data/plink/2:11723858:A_G.ped",
                   col_names = FALSE,
                   col_select = c(2, 7, 8))

colnames(snp1) <- c("ID", "a1", "a2")

snp1 <- snp1 %>% 
    mutate(snp_2_11723858_A_G  = paste0(a1, a2)) %>% 
    select(-a1, -a2)

# snp_2 - 11:30382899:C_T   cg21393004
snp2 <- read_delim("../data/plink/11:30382899:C_T.ped",
                   col_names = FALSE,
                   col_select = c(2, 7, 8))

 colnames(snp2) <- c("ID", "a1", "a2")
 
 snp2 <- snp2 %>% 
     mutate(snp_11_30382899_C_T = paste0(a1, a2)) %>% 
     select(-a1, -a2)

# snp_3 - 6:151767907:A_G   cg23088022
 
 snp3 <- read_delim("../data/plink/6:151767907:A_G.ped",
                    col_names = FALSE,
                    col_select = c(2, 7, 8))
 
 colnames(snp3) <- c("ID", "a1", "a2")
 
 snp3 <- snp3 %>% 
   mutate(snp_6_151767907_A_G = paste0(a1, a2)) %>% 
   select(-a1, -a2) 
 
 # snp_4 - 12:54387947:A_G  cg09605287
 
 snp4 <- read_delim("../data/plink/12:54387947:A_G.ped",
                    col_names = FALSE,
                    col_select = c(2, 7, 8))
 
 colnames(snp4) <- c("ID", "a1", "a2")
 
 snp4 <- snp4 %>% 
   mutate(snp_6_151767907_A_G = paste0(a1, a2)) %>% 
   select(-a1, -a2) 
 
 
# Samples to keep
d <- fread("../../../tensorQTL_final/final_20250210/data/covs_130_tensor.txt")
samples_to_keep <- colnames(d)[-1]


# Methylation data
met <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_qced.bed.gz") %>% 
  filter(CpG %in% c("cg07314298"))

met <- met %>% 
  t() %>% 
  as_tibble(rownames = "ID") %>% 
  rename_all(~ c("ID", "cg07314298")) %>% 
  filter(ID != "ID") %>%
  filter(ID %in% samples_to_keep)

data <- met %>% 
  inner_join(snp1, join_by(ID)) %>%
  mutate(cg07314298 = as.numeric(cg07314298))

# Methylation data snp_1 - 2:11723858:A_G    cg07314298
met1 <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_qced.bed.gz") %>% 
  filter(CpG %in% c("cg07314298"))

met1 <- met1 %>% 
  t() %>% 
  as_tibble(rownames = "ID") %>% 
  rename_all(~ c("ID", "cg07314298")) %>% 
  filter(ID != "ID") %>%
  filter(ID %in% samples_to_keep)

data1 <- met1 %>% 
  inner_join(snp1, join_by(ID)) %>%
  mutate(cg07314298 = as.numeric(cg07314298))

# Methylation data snp_2 - 11:30382899:C_T   cg21393004
met2 <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_qced.bed.gz") %>% 
  filter(CpG %in% c("cg21393004"))

met2 <- met2 %>% 
  t() %>% 
  as_tibble(rownames = "ID") %>% 
  rename_all(~ c("ID", "cg21393004")) %>% 
  filter(ID != "ID") %>%
  filter(ID %in% samples_to_keep)

data2 <- met2 %>% 
  inner_join(snp2, join_by(ID)) %>%
  mutate(cg21393004 = as.numeric(cg21393004))

# Methylation data snp_3 - cg23088022
met3 <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_qced.bed.gz") %>% 
  filter(CpG %in% c("cg23088022"))

met3 <- met3 %>% 
  t() %>% 
  as_tibble(rownames = "ID") %>% 
  rename_all(~ c("ID", "cg23088022")) %>% 
  filter(ID != "ID") %>%
  filter(ID %in% samples_to_keep)

data3 <- met3 %>% 
  inner_join(snp3, join_by(ID)) %>%
  mutate(cg23088022 = as.numeric(cg23088022))

# Methylation data snp_4 - cg09605287
met4 <- read_tsv("/data/projects/endometriosis/meth_QC_oct24/new_2/results/betas_qced.bed.gz") %>% 
  filter(CpG %in% c("cg09605287"))

met4 <- met4 %>% 
  t() %>% 
  as_tibble(rownames = "ID") %>% 
  rename_all(~ c("ID", "cg09605287")) %>% 
  filter(ID != "ID") %>%
  filter(ID %in% samples_to_keep)

data4 <- met4 %>% 
  inner_join(snp4, join_by(ID)) %>%
  mutate(cg09605287 = as.numeric(cg09605287))


# Methylation across genotypes ================================================
theme_std <- theme_set(theme_minimal(base_size = 11))
theme_update(
    axis.title = element_text(size = 10),
    axis.text = element_text(color = "black", size = 9),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.4),
    axis.line = element_line(linewidth = 0.4),
    axis.line.x.bottom = element_line(linewidth = 0.4),
    panel.background = element_rect(fill = 'white', color = 'white'),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(size = 10, margin = margin(b = 15)),
    # plot.tag = element_text(size = 16),
    # legend.title = element_text(size = 9),
    # legend.text = element_text(size = 8),
    # legend.key.size = unit(0.4, "cm")
)

################################################################################

# Methylation data 2:11723858:A_G    cg07314298

nums1 <- table(data1$snp_2_11723858_A_G)
nums1 <- paste0(names(nums1), "\n","(n = ", nums1, ")")

box1 <- ggplot(data1, aes(snp_2_11723858_A_G, cg07314298, fill = snp_2_11723858_A_G)) +
    geom_boxplot(width = 0.55, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "gray") +
    scale_x_discrete(labels = nums1) +
    scale_y_continuous(limits = c(0, 0.5),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c("#F4D06F", "#9DD9D2", "#FFF8F0")) +
    xlab("") +
    ylab("DNAm beta values") +
    labs(title = "cg07314298") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
    
# Export to pdf
grDevices::cairo_pdf("cg07314298 - snp_2_11723858_A_G.pdf",
                     width = 2.5, height = 4)
box1
dev.off()

# Methylation data 11:30382899:C_T   cg21393004

nums2 <- table(data$snp_2_11731846_C_T)
nums2 <- paste0(names(nums2), "\n","(n = ", nums2, ")")

box2 <- ggplot(data, aes(snp_2_11731846_C_T, cg12568421)) +
    geom_boxplot(fill = "lightgrey", width = 0.55, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    scale_x_discrete(labels = nums2) +
    scale_y_continuous(limits = c(0.2, 1),
                       expand = c(0, 0),
                       breaks = c(seq(0.2, 1, 0.2))) +
    xlab("") +
    ylab("Methylation (Beta)") +
    labs(title = "cg12568421 - rs2884374") +
    theme(legend.position = "none")


grDevices::cairo_pdf("../results/plots/boxplot_mqtl_cg12568421.pdf",
                     width = 2.5, height = 4)
box2
dev.off()


# Correlation between the two cpgs ============================================

metcorr <- cor.test(data$cg07314298, data$cg12568421)
metcorr1 <- signif(metcorr$estimate, 3)
metcorr2 <- signif(metcorr$p.value, 3)

textbox <- paste0(
    "Pearson corr. estimate = ", metcorr1, "\n",
    "P-value = ", metcorr2, "\n"
)


c1 <- ggplot(data, aes(cg07314298, cg12568421)) +
    geom_smooth(method = lm, se = F, alpha = 0.1, color = "tomato3") +
    geom_point(size = 1) +
    scale_x_continuous(limits = c(0, 0.5),
                       expand = c(0, 0),
                       breaks = seq(0, 0.5, 0.1)) +
    scale_y_continuous(limits = c(0.3, 0.9),
                       expand = c(0, 0),
                       breaks = seq(0.3, 0.9, 0.1)) +
    xlab("cg07314298") +
    ylab("cg12568421") +
    # labs(title = "") +
    annotate(geom = "text", label = textbox, x = 0.35, y = 0.5,
             size = 3) +
    theme(legend.position = "none",
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.line.x.bottom = element_line(linewidth = 0.4),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))

grDevices::cairo_pdf("../results/plots/corr_cg07314298_cg12568421.pdf",
                     width = 4, height = 4)
c1
dev.off()


# Methylation levels 14 controls vs 14 cases ==================================

# Read previously generated methylome file
metbed <- read_tsv("../../mqtl_efa_01/tensorQTL/textfiles/methylome_BED.txt")

# Remove controls 75 samples provided by Iraia,
# keeping 14 cases and 15 controls
samples <- read_lines("../data/samples_discovery_mqtl.txt")

metbed <- metbed %>% 
    select(!all_of(samples))

# Methylation data
met2 <- metbed %>% 
    filter(ID %in% c("cg07314298", "cg12568421")) %>% 
    select(-`#Chr`, -start, -end)

met2 <- met2 %>% 
    t() %>% 
    as_tibble(rownames = "ID") %>% 
    rename_all(~ c("ID", "cg07314298", "cg12568421")) %>% 
    slice(-1)

met2 <- met2 %>% 
    mutate(cg07314298 = as.numeric(cg07314298)) %>% 
    mutate(cg12568421 = as.numeric(cg12568421)) %>% 
    mutate(Group = ifelse(str_detect(ID, "ENDOCONT"), "Control", "Case")) %>% 
    mutate(Group = factor(Group, levels = c("Control", "Case")))

numsvs <- table(met2$Group)
numsvs <- paste0(names(numsvs), "\n","(n = ", numsvs, ")")

boxvs1 <- ggplot(met2, aes(Group, cg07314298)) +
    geom_boxplot(fill = "lightgrey", width = 0.5, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    scale_x_discrete(labels = numsvs) +
    scale_y_continuous(limits = c(0, 0.4),
                       expand = c(0, 0),
                       breaks = c(seq(0, 0.4, 0.1))) +
    xlab("") +
    ylab("Methylation (Beta)") +
    labs(title = "cg07314298") +
    theme(legend.position = "none")

grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg07314298.pdf",
                     width = 2, height = 3.5)
boxvs1
dev.off()

boxvs2 <- ggplot(met2, aes(Group, cg12568421)) +
    geom_boxplot(fill = "lightgrey", width = 0.5, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    scale_x_discrete(labels = numsvs) +
    scale_y_continuous(limits = c(0.5, 0.9),
                       expand = c(0, 0),
                       breaks = c(seq(0.5, 0.9, 0.1))) +
    xlab("") +
    ylab("Methylation (Beta)") +
    labs(title = "cg12568421") +
    theme(legend.position = "none")


grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg12568421.pdf",
                     width = 2, height = 3.5)
boxvs2
dev.off()


# Also add genotype data ======================================================

# system("plink --bfile /data/projects/endometriosis/imputation/imputed-rsq09/whole_genome_maf05_hw0.05 --snp 2:11723110:G:A --recode --out ../data/plink2/2_11723110_A_G")
# system("plink --bfile /data/projects/endometriosis/imputation/imputed-rsq09/whole_genome_maf05_hw0.05 --snp 2:11731846:T:C --recode --out ../data/plink2/2_11731846_C_T")

snp1 <- read_delim("../data/plink2/2_11723110_A_G.ped",
                   col_names = FALSE,
                   col_select = c(2, 7, 8))

colnames(snp1) <- c("ID", "a1", "a2")

snp1 <- snp1 %>% 
    mutate(snp_2_11723110_A_G = paste0(a1, a2)) %>% 
    select(-a1, -a2)

snp2 <- read_delim("../data/plink2/2_11731846_C_T.ped",
                   col_names = FALSE,
                   col_select = c(2, 7, 8))

colnames(snp2) <- c("ID", "a1", "a2")

snp2 <- snp2 %>% 
    mutate(snp_2_11731846_C_T = paste0(a1, a2)) %>% 
    select(-a1, -a2)


met2 <- met2 %>% 
    inner_join(snp1, join_by(ID)) %>% 
    inner_join(snp2, join_by(ID))


boxgm1 <- ggplot(met2, aes(Group, cg07314298, fill = snp_2_11723110_A_G)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0, position = position_dodge(0.8)) +
    geom_jitter(size = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
    scale_x_discrete(labels = numsvs) +
    scale_y_continuous(limits = c(0, 0.4),
                       expand = c(0, 0),
                       breaks = c(seq(0, 0.4, 0.1))) +
    rcartocolor::scale_fill_carto_d(palette = "Pastel", name = "") +
    xlab("") +
    ylab("Methylation (Beta)") +
    labs(title = "cg07314298 - rs16857668") +
    theme(legend.position = "bottom")

grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg07314298_snp_2_11723110_A_G.pdf",
                     width = 4, height = 5)
boxgm1
dev.off()


boxgm2 <- ggplot(met2, aes(Group, cg12568421, fill = snp_2_11731846_C_T)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0, position = position_dodge(0.8)) +
    geom_jitter(size = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
    scale_x_discrete(labels = numsvs) +
    scale_y_continuous(limits = c(0.5, 0.9),
                       expand = c(0, 0),
                       breaks = c(seq(0.5, 0.9, 0.1))) +
    rcartocolor::scale_fill_carto_d(palette = "Pastel", name = "") +
    xlab("") +
    ylab("Methylation (Beta)") +
    labs(title = "cg12568421 - rs2884374") +
    theme(legend.position = "bottom")

grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg12568421_snp_2_11731846_C_T.pdf",
                     width = 4, height = 5)
boxgm2
dev.off()


# Alternative version
# 
# boxgm1alt <- ggplot(met2, aes(Group, cg07314298, color = snp_2_11723110_A_G)) +
#     geom_jitter(size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
#     scale_x_discrete(labels = numsvs) +
#     scale_y_continuous(limits = c(0, 0.4),
#                        expand = c(0, 0),
#                        breaks = c(seq(0, 0.4, 0.1))) +
#     rcartocolor::scale_color_carto_d(palette = "Vivid", name = "") +
#     xlab("") +
#     ylab("Methylation (Beta)") +
#     labs(title = "cg07314298 - rs16857668") +
#     theme(legend.position = "bottom")
# 
# grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg07314298_snp_2_11723110_A_G_v2.pdf",
#                      width = 4, height = 5)
# boxgm1alt
# dev.off()
# 
# 
# boxgm2alt <- ggplot(met2, aes(Group, cg12568421, color = snp_2_11731846_C_T)) +
#     geom_jitter(size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
#     scale_x_discrete(labels = numsvs) +
#     scale_y_continuous(limits = c(0.5, 0.9),
#                        expand = c(0, 0),
#                        breaks = c(seq(0.5, 0.9, 0.1))) +
#     rcartocolor::scale_color_carto_d(palette = "Vivid", name = "") +
#     xlab("") +
#     ylab("Methylation (Beta)") +
#     labs(title = "cg12568421 - rs2884374") +
#     theme(legend.position = "bottom")
# 
# grDevices::cairo_pdf("../results/plots/boxplot_caseVScont_cg12568421_snp_2_11731846_C_T_v2.pdf",
#                      width = 4, height = 5)
# boxgm2alt
# dev.off()


# GREB1 expression ============================================================

greb <- read_tsv("../data/expression_greb1.csv")

greb <- greb %>% 
    rename("rel_exp" = `Relative expression GREB1`) %>% 
    mutate(rel_exp = str_replace_all(rel_exp, ",", ".")) %>% 
    mutate(rel_exp = as.numeric(rel_exp)) %>% 
    filter(Sample != "ENDOCONT034") %>% 
    filter(!is.na(rel_exp))

met3 <- met2 %>% 
    inner_join(greb, join_by(ID == Sample))


corcg1 <- cor.test(met3$cg07314298, met3$rel_exp)
corcg1_e <- signif(corcg1$estimate, 3)
corcg1_p <- signif(corcg1$p.value, 3)

textbox_cg1 <- paste0(
    "Pearson corr. estimate = ", corcg1_e, "\n",
    "P-value = ", corcg1_p, "\n"
)

exp_cg1 <- ggplot(met3, aes(cg07314298, rel_exp)) +
    geom_smooth(method = lm, se = F, alpha = 0.1, color = "tomato3") +
    geom_point(size = 1) +
    scale_x_continuous(limits = c(0.1, 0.4),
                       expand = c(0, 0),
                       breaks = seq(0.1, 0.5, 0.1)) +
    scale_y_continuous(limits = c(0.9, 1.6),
                       expand = c(0, 0),
                       breaks = seq(0.9, 1.6, 0.1)) +
    xlab("cg07314298") +
    ylab("GREB1 Relative Expression") +
    labs(title = "") +
    annotate(geom = "text", label = textbox_cg1, x = 0.2, y = 1.3, size = 3) +
    theme(legend.position = "none",
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.line.x.bottom = element_line(linewidth = 0.4),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))

grDevices::cairo_pdf("../results/plots/corr_greb1_cg07314298.pdf",
                     width = 4, height = 4)
exp_cg1
dev.off()


corcg2 <- cor.test(met3$cg12568421, met3$rel_exp)
corcg2_e <- signif(corcg2$estimate, 3)
corcg2_p <- signif(corcg2$p.value, 3)

textbox_cg2 <- paste0(
    "Pearson corr. estimate = ", corcg2_e, "\n",
    "P-value = ", corcg2_p, "\n"
)

exp_cg2 <- ggplot(met3, aes(cg12568421, rel_exp)) +
    geom_smooth(method = lm, se = F, alpha = 0.1, color = "tomato3") +
    geom_point(size = 1) +
    scale_x_continuous(limits = c(0.6, 0.9),
                       expand = c(0, 0),
                       breaks = seq(0.5, 0.9, 0.1)) +
    scale_y_continuous(limits = c(0.9, 1.6),
                       expand = c(0, 0),
                       breaks = seq(0.9, 1.6, 0.1)) +
    xlab("cg12568421") +
    ylab("GREB1 Relative Expression") +
    labs(title = "") +
    annotate(geom = "text", label = textbox_cg2, x = 0.70, y = 1.25, size = 3) +
    theme(legend.position = "none",
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.line.x.bottom = element_line(linewidth = 0.4),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))

grDevices::cairo_pdf("../results/plots/corr_greb1_cg12568421.pdf",
                     width = 4, height = 4)
exp_cg2
dev.off()


# GREB1 and SNP genotype

met3 <- met3 %>% 
    inner_join(snp1, join_by(ID)) %>% 
    inner_join(snp2, join_by(ID))
    
nums1 <- table(met3$snp_2_11723110_A_G)
nums1 <- paste0(names(nums1), "\n","(n = ", nums1, ")")

nums2 <- table(met3$snp_2_11731846_C_T)
nums2 <- paste0(names(nums2), "\n","(n = ", nums2, ")")

exp_snp1 <- ggplot(met3, aes(snp_2_11723110_A_G, rel_exp)) +
    geom_boxplot(fill = "lightgrey", width = 0.55, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    scale_x_discrete(labels = nums1) +
    scale_y_continuous(limits = c(0.9, 1.6),
                       expand = c(0, 0),
                       breaks = seq(0.9, 1.6, 0.1)) +
    xlab("rs16857668") +
    ylab("GREB1 Relative Expression") +
    theme(legend.position = "none")


grDevices::cairo_pdf("../results/plots/boxplot_greb1_snp_2_11723110_A_G.pdf",
                     width = 2, height = 3.5)
exp_snp1
dev.off()


exp_snp2 <- ggplot(met3, aes(snp_2_11731846_C_T, rel_exp)) +
    geom_boxplot(fill = "lightgrey", width = 0.55, outlier.alpha = 0) +
    geom_jitter(width = 0.15, fill = "grey", size = 0.8) +
    scale_x_discrete(labels = nums2) +
    scale_y_continuous(limits = c(0.9, 1.6),
                       expand = c(0, 0),
                       breaks = seq(0.9, 1.6, 0.1)) +
    xlab("") +
    xlab("rs2884374") +
    ylab("GREB1 Relative Expression") +
    theme(legend.position = "none")


grDevices::cairo_pdf("../results/plots/boxplot_greb1_snp_2_11731846_C_T.pdf",
                     width = 2.5, height = 3.5)
exp_snp2
dev.off()


# GREB1 - Case control (Germany)

greb_g <- greb %>% 
    filter(Source != "EFA") %>% 
    mutate(Condition = factor(Condition, levels = c("Control", "Case")))
    
nums_g <- table(greb_g$Condition)
nums_g <- paste0(names(nums_g), "\n","(n = ", nums_g, ")")


greb_case <- ggplot(greb_g, aes(Condition, rel_exp)) +
    geom_boxplot(fill = "lightgrey", width = 0.55, outlier.alpha = 0) +
    scale_x_discrete(labels = nums_g) +
    geom_jitter(width = 0.15, size = 0.8) +
    scale_y_continuous(limits = c(0.9, 1.6),
                       expand = c(0, 0),
                       breaks = seq(0.9, 1.6, 0.1)) +
    xlab("") +
    ylab("GREB1 Relative Expression") +
    theme(legend.position = "none")


grDevices::cairo_pdf("../results/plots/boxplot_greb1_controlVScase.pdf",
                     width = 2, height = 3.5)
greb_case
dev.off()

