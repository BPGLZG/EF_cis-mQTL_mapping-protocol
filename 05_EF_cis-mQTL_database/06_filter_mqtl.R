# February 10, 2025
#06_filter_mqtl.R

#Once tensor has been executed, concatenate all parquets and filter by <5e-08

# Set working directory
setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210")

#!/usr/bin/env Rscript

#Sys.setenv("NOT_CRAN" = "true", "LIBARROW_BUILD" = FALSE, "ARROW_R_DEV" = TRUE)
#install.packages("arrow")

library(tidyverse)
library(data.table)
library(arrow)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load EPIC manifest
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>% 
  as_tibble() %>% 
  select(Name, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group)

# Read parquet, filter by p-value and add the Gene annotation from EPIC
# manifest for a quicker gene identification
filt_mqtl <- function(filename, met_anno) {
  
  data <- filename %>% 
    read_parquet() %>% 
    filter(pval_nominal < 5e-08) %>% 
    inner_join(met_anno, by = join_by(phenotype_id == Name))
  
  return(data)
  
}

# List of parquet files
filename_list <- dir(path = "/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_nom",
                     pattern = "\\.parquet", 
                     full.names = TRUE)

# Apply the function to each element of the list
filtered_mqtl <- map(filename_list, filt_mqtl, met_anno = anno, .progress = TRUE)

# Bind the results into a single tibble
filtered_mqtl <- filtered_mqtl %>% 
  bind_rows() %>% 
  arrange(variant_id)

# Write
write_tsv(filtered_mqtl, "/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_nom/tensorqtl_130_5e-08_20250210.tsv")

# Load results
tensor <- fread("/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_nom/tensorqtl_130_5e-08_20250210.tsv")
dim(tensor)   # 1190678 12
tensor[1:5,1:12]



alleles <- map2(formatted[, EA], formatted[, OA], sort_and_paste)
formatted[, CPAID := paste0(CHR, ":", POS, ":", alleles)]
[13:10, 21/5/2024] Pato Gonzalez: 