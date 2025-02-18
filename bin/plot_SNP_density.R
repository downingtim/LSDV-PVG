#!/usr/bin/env Rscript
library(ggplot2)
options(warn=-1)  # Suppress warnings
x <- read.csv("gfavariants.vcf", sep = "\t", skip = 4)
x$POS <- x$POS / 1000

hist.data <- hist(x$POS, plot = FALSE, breaks = 9999)
hist.df <- data.frame(
  mids = hist.data$mids,
  counts = log10(hist.data$counts + 1) )

p <- ggplot(hist.df, aes(x = mids, y = counts)) +
  geom_bar(stat = "identity", fill = "grey", alpha=0.6) +
  xlab("Genomic position (Kb)") +
  ylab("Log10(Mutation count + 1)") +
  theme_minimal()

# Save the plot to a PDF file
ggsave("mutation_density.pdf", plot = p, width = 42, height = 6)