#!/usr/bin/env Rscript

# make heatmap
x2 <- as.data.frame(read.csv("../../../CURRENT/PANAROO/gene_presence_absence.Rtab", sep="\t"))
# str(x2) # samples as columns, genes as rows
rownames(x2) <- x2[,1] # reassign row (ie, gene) names
x2 <- x2[,-1] # remove gene names from table (not needed)
my_palette <- colorRampPalette(c("lightgray", "green", "black"))(n = 299)

library(gplots)

# x2[x2 == 0] <- NA

pdf("../../../CURRENT/PANAROO/heatmap.pdf", width=20, height=40)
col_breaks = c(seq(0.001,0.2,length=100), seq(0.3,0.4,length=100), seq(0.5,0.6,length=100))
heatmap.2(as.matrix(x2), dendrogram = "none", col=my_palette, margins = c(59,8),
trace = "none", density.info = "none", breaks =  col_breaks, na.color = "white")
dev.off()

q()
y
