# imiss-vs-het.R script

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("geneplotter")

setwd('/data/projects/endometriosis/qc_final')

imiss = read.table('all.imiss', header = T)
imiss$logF_MISS = log10(imiss[,6])

het = read.table('plink.het', header = T)
het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.

colors <- densCols(imiss$logF_MISS, het$meanHet)
pdf('imiss-vs-het.pdf')

plot(imiss$logF_MISS, het$meanHet, col = colors, xlim = c(-3,0), ylim = c(0.0,5), pch = 20, xlab = 'Proportion of missing genotypes', ylab = 'Heterozygosity rate', axes = F)
axis(2, at = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50), tick = T)
axis(1, at = c(-3, -2, -1, 0), labels = c(0.001, 0.01, 0.1, 1))
abline(h = mean(het$meanHet) - (4*sd(het$meanHet)), col = 'RED', lty =2)
abline(h = mean(het$meanHet) + (4*sd(het$meanHet)), col = 'RED', lty =2)
abline(v = -1.522879, col = 'RED', lty = 2)

dev.off()