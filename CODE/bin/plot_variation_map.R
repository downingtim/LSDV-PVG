#!/usr/bin/env Rscript

library(ggplot2)
library(ggExtra)
map <- read.csv("variation_map.txt", sep="\t")
colnames(map) <- c("sample", "allele1", "pos1", "allele2", "pos2")
pdf("variation_map-basic.pdf")
plot(map$pos1, map$pos2, cex=0.05) 
grid()
dev.off()
mapd <- subset(map, abs(map$pos1-map$pos2) >600)
unique(mapd$sample)

p1 <- ggplot(map, aes(x =pos1, y = log10( abs(pos1-pos2)), fill = sample, color=sample)) +
geom_point(size=0.1, show.legend = F) + 
geom_line(size=0.6, alpha=0.8, linetype=4, show.legend = F) +
labs(y= "Log10-scaled other sample genome position", x="Main sample genome position") +
theme(legend.position = "none") 
pdf("variation_map.pdf", width=13, height=6)
ggExtra::ggMarginal(p1, type="densigram", size=5, colour="grey", bins=50)  
dev.off()

