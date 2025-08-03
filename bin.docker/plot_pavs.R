#!/usr/bin/env Rscript

# Load required libraries with error handling
required_packages <- c("data.table", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <input_file.tsv> [output_file.pdf]\n")
  cat("Example: Rscript script.R data.tsv output.pdf\n")
  quit(status = 1)
}

input_file <- as.character(args[1])
output_file <- if (length(args) >= 2) as.character(args[2]) else gsub("tsv", "pdf", as.character(args[1]))

# Validate input file exists
if (!file.exists(input_file)) {
  cat("Error: Input file", input_file, "does not exist.\n")
  quit(status = 1)
}

cat("Reading data from:", input_file, "\n")
start_time <- Sys.time()

# Fast data reading with data.table
x <- fread(input_file, sep = "\t", select = "pav", showProgress = FALSE)

# Check if 'pav' column exists
if (!"pav" %in% names(x)) {
  cat("Error: Column 'pav' not found in the input file.\n")
  quit(status = 1)
}

# Remove NA values efficiently
x <- x[!is.na(pav)]
n_samples <- nrow(x)

cat("Data loaded successfully. Samples:", n_samples, "\n")
cat("Data range:", min(x$pav), "to", max(x$pav), "\n")

# Calculate optimal number of bins for large datasets
optimal_breaks <- min(62, max(30, ceiling(sqrt(n_samples))))
cat("Using", optimal_breaks, "bins for histogram\n")

# Create enhanced plot
cat("Generating plot...\n")
pdf(output_file, width = 7, height = 6)

# Set up plotting parameters for better appearance
par(mar = c(5, 5, 4, 2) + 0.1,  # Increase margins
    bg = "white",
    family = "sans")

# Create histogram with enhanced styling
h <- hist(x$pav, 
          breaks = optimal_breaks,
          xlab = "Number of samples",
          ylab = "Frequency",
          main = paste("Distribution of PAV values (n =", format(n_samples, big.mark = ","), "PAVs)"),
          freq = FALSE,
          col = "#4292C6",           # Professional blue color
          border = "#08519C",        # Darker blue border
          las = 1,                   # Horizontal y-axis labels
          cex.lab = 1.2,            # Larger axis labels
          cex.main = 1.3,           # Larger title
          cex.axis = 1.1)           # Larger axis numbers

# Add grid for better readability
grid(col = "lightgray", lty = "dotted", lwd = 1)

# Add histogram back on top of grid
hist(x$pav, 
     breaks = optimal_breaks,
     freq = FALSE,
     col = "#4292C6",
     border = "#08519C",
     add = TRUE)

# Add density curve overlay
if (n_samples > 1000) {  # Only for reasonably large datasets
  tryCatch({
    # Sample data for density estimation if very large
    sample_size <- min(50000, n_samples)
    sample_data <- sample(x$pav, sample_size)
    dens <- density(sample_data, na.rm = TRUE)
    lines(dens, col = "#E31A1C", lwd = 2)
  }, error = function(e) {
    cat("Note: Could not add density curve\n")
  })
}

# Add summary statistics as text
mean_val <- mean(x$pav, na.rm = TRUE)
median_val <- median(x$pav, na.rm = TRUE)
sd_val <- sd(x$pav, na.rm = TRUE)

# Position text in upper left
legend("topleft", 
       legend = c(paste("Mean:", round(mean_val, 2)),
                  paste("Median:", round(median_val, 2)),
                  paste("SD:", round(sd_val, 2))),
       bty = "n",  # No box around legend
       cex = 1.3,
       text.col = "#333333")

dev.off()

end_time <- Sys.time()
execution_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

cat("Plot saved to:", output_file, "\n")
cat("Execution time:", execution_time, "seconds\n")
cat("Memory usage:", format(object.size(x), units = "MB"), "\n")

