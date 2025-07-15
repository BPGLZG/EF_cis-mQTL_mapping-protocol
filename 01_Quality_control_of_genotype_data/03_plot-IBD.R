# plot-IBD.R script

setwd("/data/projects/endometriosis/qc_final")

data = read.table("all_geno0.05_hwe1e06_maf0.005_ind_IBD.genome", h = T)

pdf("all_geno0.05_hwe1e06_maf0.005_ind_IBD_hist.pdf")
hist(data$PI_HAT, ylim = c(0,100), col = "RED", breaks = 100, xlab = "Estimated mean pairwise IBD", main = "")

dev.off()

out<-data[data$PI_HAT>0.185,]

write.table(out,"all_geno0.05_hwe1e06_maf0.005_ind_IBD_fail_IBD_check.txt", sep = "\t", quote = F, row.names = F)

quit ()  
