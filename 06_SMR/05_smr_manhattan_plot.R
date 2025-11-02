##### Manhattan plot - SMR hits
### February 27, 2025

### Manhattan plot 1) GCST90205183  - Endometriosis
### Manhattan plot 2) GCST90269970  - Endometriosis
### Manhattan plot 3) GCST010098    - Breast cancer
### Manhattan plot 4) GST006464     - Endometrial cancer
### Manhattan plot 5) GCST90244167  - Ovarian cancer

### Set working directory
setwd("/data/projects/endometriosis/smr_final/smr_130_5_9/plots/")

### Install packages
install.packages('qqman')

# Load libraries
library('qqman')
library(data.table)

################################################################################
### Manhattan plot 1) GCST90205183 - Endometriosis
# Load results
endo <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90205183/GCST90205183.msmr")

endo_heidi <- endo[p_HEIDI > 0.05,]

# Plot GCST90205183
png("GCST90205183_manhattan_plot.png",
    width = 10,
    height = 6,
    units = "in",
    res = 300
)

manhattan(endo_heidi,
          chr = 'ProbeChr',
          bp = 'Probe_bp',
          p = 'p_SMR_multi',
          snp = 'probeID',
          col = c('gray10', 'gray60'),
          genomewideline = -log10(1.992905e-06),
          suggestiveline = FALSE
          
    )

dev.off()

################################################################################
### Manhattan plot 2) GCST90269970 - Endometriosis
# Load results
endo_2 <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometriosis/GCST90269970/GCST90269970.msmr")

endo_2_heidi <- endo_2[p_HEIDI > 0.05,]

endo_2_heidi <- endo_2_heidi[!is.na(p_SMR_multi), ]

# Plot GCST90269970
png("GCST90269970_manhattan_plot.png",
    width = 10,
    height = 6,
    units = "in",
    res = 300
)

manhattan(endo_2_heidi,
          chr = "ProbeChr",
          bp = "Probe_bp",
          p = "p_SMR_multi",
          snp = "probeID",
          col = c("gray10", "gray60"),
          genomewideline = -log10(1.992905e-06),
          suggestiveline = FALSE
)

dev.off()

################################################################################
### Manhattan plot 3) GCST010098 - Breast cancer
# Load results
breast <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/breast_cancer/GCST010098/GCST010098.msmr")

breast_heidi <- breast[p_HEIDI > 0.05,]

breast_heidi <- breast_heidi[!is.na(p_SMR_multi), ]

# Plot GCST010098
png('GCST010098_manhattan_plot.png',
    width = 10,
    height = 6,
    units = "in",
    res = 300
)

manhattan(breast_heidi,
          chr = 'ProbeChr',
          bp = 'Probe_bp',
          p = 'p_SMR_multi',
          snp ='probeID',
          col = c('gray10', 'gray60'),
          genomewideline = -log10(1.992905e-06),
          suggestiveline = FALSE
)

dev.off()

################################################################################
### Manhattan plot 4) GST006464 - Endometrial cancer
# Load results
endometrial <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/endometrial_cancer/GCST006464/GCST006464.msmr")

endometrial_heidi <- endometrial[p_HEIDI > 0.05,]

endometrial_heidi <- endometrial_heidi[!is.na(p_SMR_multi), ]

# Plot GST006464
png('GCST006464_manhattan_plot.png',
    width = 10,
    height = 6,
    units = "in",
    res = 300
)

manhattan(endometrial_heidi,
          chr = 'ProbeChr',
          bp = 'Probe_bp',
          p = 'p_SMR_multi',
          snp ='probeID',
          col = c('gray10', 'gray60'),
          genomewideline = -log10(1.992905e-06),
          suggestiveline = FALSE
)

dev.off()

################################################################################
### Manhattan plot 5) GCST90244167  - Ovarian cancer
# Load results
ovarian <- fread("/data/projects/endometriosis/smr_final/smr_130_5_9/results/ovarian_cancer/ieu-a-1120/ieu-a-1120.msmr")

ovarian_heidi <- ovarian[p_HEIDI > 0.05,]

ovarian_heidi <- ovarian_heidi[!is.na(p_SMR_multi), ]

# Plot GCST90244167
png('GCST90244167_manhattan_plot.png',
    width = 10,
    height = 6,
    units = "in",
    res = 300
)

manhattan(ovarian_heidi,
          chr = 'ProbeChr',
          bp = 'Probe_bp',
          p = 'p_SMR_multi',
          snp ='probeID',
          col = c('gray10', 'gray60'),
          genomewideline = -log10(1.992905e-06),
          suggestiveline = FALSE
)

dev.off()
