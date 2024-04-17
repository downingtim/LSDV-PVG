#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
name=args[1]
x <- read.csv(name, sep="\t", skip=4)
#str(x)
x$POS <- x$POS/1000
hist.data = hist(x$POS, plot=F, breaks=9999)
hist.data$counts = log10(hist.data$counts  + 1) 
pdf("../../../CURRENT/mutation_density.pdf", width=20, height=8)
plot( hist.data, xlab="Genomic position (Kb)", col="grey")
dev.off()