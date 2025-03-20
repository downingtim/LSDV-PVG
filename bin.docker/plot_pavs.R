#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <arg1> \n")
  quit(status = 1) }

x <- read.csv(args[1], sep="\t")
new <- gsub("tsv", "pdf", args[1])

pdf(new, width=8, height=7)
hist( x$pav, breaks=62, xlab=c("Number of samples"),  ylab=c("Number of mutations"), main="", freq=F)  
dev.off()
