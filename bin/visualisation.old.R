#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <arg1> \n")
  quit(status = 1) }

arg <- as.numeric(args[1])
library(ggplot2)
h1 <- read.csv("heaps.txt", sep="\t", header=F) 
colnames(h1) <- c("a", "b", "pangenome")
valu = 2/arg
str(h1)

pdf("heaps.pdf", width=11, height=8)
ggplot()  + geom_point(aes(x=arg:valu,y=h1$pangenome), color="red", size=1.5, alpha=0.9)  +
	  geom_line(aes(x=arg:valu, y=h1$pangenome),  col="red",  linewidth=1.5, alpha=0.9) +
	  scale_y_continuous(limits = c(min(h1$pangenome)-1000, max(h1$pangenome)+1000))  +
	  labs(x = "Fraction of samples") +
	  labs(y = "Number of bases in paths") + theme(legend.position="top") + guides()
dev.off()
