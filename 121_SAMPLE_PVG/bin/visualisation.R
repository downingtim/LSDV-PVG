#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <arg1> \n")
  quit(status = 1) }

arg <- as.numeric(args[1])
print(arg)
library(ggplot2)

output_file <- file("../../../data/heaps.txt", "w")

for (e in 1:arg){ 
  args1 <- c(paste("heaps -i ../../../CURRENT/*.og -n 1000 -d ", e, " -t 12 | sort -nk 2 | tail -n 1", sep=""))
  output <- system2(command="odgi", args=args1, stdout = TRUE, stderr = FALSE)
  out2 <- c(output, "\n")
  cat(out2, file = output_file, append = TRUE)
}

close(output_file)

h1 <- read.csv("../../../data/heaps.txt", sep="\t", header=F) 
colnames(h1) <- c("a", "b", "pangenome")
valu = 2/arg
str(h1)

pdf("../../../CURRENT/heaps.pdf", width=11, height=8)
ggplot()  + geom_point(aes(x=arg:valu,y=h1$pangenome), color="red", size=1.5, alpha=0.9)  +
	  geom_line(aes(x=arg:valu, y=h1$pangenome),  col="red",  linewidth=1.5, alpha=0.9) +
	  scale_y_continuous(limits = c(min(h1$pangenome)-1000, max(h1$pangenome)+1000))  +
	  labs(x = "Fraction of samples") +
	  labs(y = "Number of bases in paths") + theme(legend.position="top") + guides()
dev.off()