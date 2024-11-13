#!/usr/bin/env Rscript

xx <- read.csv("lumpy_skin_disease_2_virus_2_2_.fasta.paf.txt", sep="\t")
y <- xx$num.X

pdf("../../../CURRENT/PAF_plot.pdf", width=20, height=8)
hist(y[y>0], breaks=199, col="grey", xlab="Number of mismatches between sample pairs", ylab="Number of sample pairs", main="") 
dev.off()

