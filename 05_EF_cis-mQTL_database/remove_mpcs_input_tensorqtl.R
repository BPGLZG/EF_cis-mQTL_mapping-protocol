# Set working directory
setwd("/data/projects/endometriosis/tensorQTL_final/final_20250210")

# Load required packages
library(readr)

# Remove mPCS (8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20)

# Load (20) residualized mPCs
mpcs <- read_tsv("data/residualized_mPCs.tsv", 
                 col_types = cols(.default = col_character()))
dim(mpcs)      # 130 21

# Remove mPCs
mpcs_rem <- mpcs %>% select(-mPC8, -mPC9, -mPC12, -mPC13, -mPC14, -mPC15, -mPC16, -mPC17, -mPC18, -mPC19, -mPC20)

dim(mpcs_rem)  # 130  10

# Transpose file
mpcs_t <- t(mpcs_rem) %>%
  as.data.frame()

# Save (9) residualized mPCs
write.table(mpcs_t, "data/residualized_11_mPC.transposed.tsv",
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
