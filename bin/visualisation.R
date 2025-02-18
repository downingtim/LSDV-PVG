#!/usr/bin/env Rscript
options(warn=-1)  # Suppress warnings

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <arg1>\n")
  quit(status = 1)
}

arg <- as.numeric(args[1])
library(ggplot2)
library(scales)  # Ensure scales is loaded

# Read the data
h1 <- read.csv("heaps.txt", sep="\t", header=F)
colnames(h1) <- c("a", "b", "pangenome")

# Generate fractions for x-axis labels
n_samples <- nrow(h1)
x_values <- seq(1/arg, 1, length.out = n_samples)
fraction_labels <- round(seq(1/arg, 1, 1/arg),2)

pdf("heaps.pdf", width=7, height=5)
ggplot(h1, aes(x = x_values, y = pangenome)) +
  geom_point(color="red", size=2, alpha=0.9) +
  geom_line(color="red", linewidth=2, alpha=0.9) +
  scale_x_continuous(breaks = x_values, labels = fraction_labels) +
  scale_y_continuous(limits = c(min(h1$pangenome) - 1000, max(h1$pangenome) + 1000)) +
  labs(x = "Fraction of samples", y = "Number of bases in paths") +
  theme_minimal()
dev.off()
