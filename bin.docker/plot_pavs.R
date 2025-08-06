#!/usr/bin/env Rscript

# Load required libraries with error handling
required_packages <- c("ggplot2","dplyr")
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

# Fast data reading with proper error handling
tryCatch({
  x <- read.csv(input_file, sep="\t", header = FALSE)
  colnames(x) <- c("value","count")
  
  # Ensure data types are correct
  x$value <- as.numeric(x$value)
  x$count <- as.numeric(x$count)
  
  # Remove any rows with NA values
  x <- x[complete.cases(x), ]
  
  if (nrow(x) == 0) {
    cat("Error: No valid data rows found in input file\n")
    quit(status = 1)
  }
  
}, error = function(e) {
  cat("Error reading input file:", e$message, "\n")
  quit(status = 1)
})

n_samples <- sum(x$count)

cat("Data loaded successfully. Samples:", n_samples, "\n")
cat("Data range:", min(x$value), "to", max(x$value), "\n")
cat("Number of unique values:", nrow(x), "\n")

# Calculate optimal number of bins for large datasets
optimal_breaks <- min(62, max(30, ceiling(sqrt(n_samples))))
cat("Using", optimal_breaks, "bins for histogram\n")
breaks <- seq(min(x$value), max(x$value), length.out = optimal_breaks + 1)

# Create enhanced plot
cat("Generating plot...\n")
pdf(output_file, width = 7, height = 6)

# Set up plotting parameters
par(mar = c(5, 5, 4, 2) + 0.1,  # Increase margins
    bg = "white",
    family = "sans")

# Assign counts to intervals
x$interval <- findInterval(x = x$value, vec = breaks, rightmost.closed = FALSE, left.open = FALSE)
x.freq <- x %>% group_by(interval) %>% summarize(total_count = sum(count), .groups = 'drop')
colnames(x.freq) <- c("interval","count")

# Design a combined table with all intervals
counts.filled <- data.frame(interval = 1:optimal_breaks, count = 0)
result <- counts.filled %>% 
  left_join(x.freq, by = "interval") %>% 
  select(interval, count.y)
result$count.y[is.na(result$count.y)] <- 0

# Create barplot
bp <- barplot(result$count.y,
              names.arg = round(breaks[-length(breaks)], 3),  # Remove last break for labeling
              xlab = "PAV values",
              ylab = "Number of mutations",
              main = paste("Distribution of PAV values (n =", format(n_samples, big.mark = ","), "PAVs)"),
              col = "#4292C6",
              border = "#08519C",
              las = 1,
              cex.lab = 1.2,
              cex.main = 1.3,
              cex.axis = 0.8)  # Smaller axis text to prevent overlap

# Add grid
grid(col = "lightgray", lty = "dotted", lwd = 1)

# Get the plotting region coordinates
usr <- par("usr")
plot_xlim <- usr[1:2]

# Generate density curve for large datasets
if (n_samples > 1000) {  
  tryCatch({
    # Sample data for density estimation if very large
    sample_size <- min(50000, n_samples)
    probabilities <- x$count / sum(x$count)
    sample_data <- sample(x$value,
                        size = sample_size,
                        replace = TRUE,
                        prob = probabilities)
    dens <- density(sample_data, na.rm = TRUE)

    # Scale y to match histogram scale
    scale_factor <- max(result$count.y) / max(dens$y)
    dens$y <- dens$y * scale_factor

    # Scale x to plot coordinates
    data_range <- range(breaks)
    dens_x_scaled <- (dens$x - data_range[1]) / diff(data_range) * diff(plot_xlim) + plot_xlim[1]

    # Add density line
    lines(dens_x_scaled, dens$y, col = "#E31A1C", lwd = 2)

  }, error = function(e) {
    cat("Note: Could not add density curve:", e$message, "\n")
  })
}

# Calculate summary statistics
mean_val <- sum(x$value * x$count) / sum(x$count)

# Calculate median
x_sorted <- x[order(x$value), ]
cumulative_count <- cumsum(x_sorted$count)
half_total <- sum(x$count) / 2
median_index <- which(cumulative_count >= half_total)[1]
median_val <- x_sorted$value[median_index]

# Calculate standard deviation
variance <- sum(x$count * (x$value - mean_val)^2) / sum(x$count)
sd_val <- sqrt(variance)

# Position text in upper left
legend("topleft",
       legend = c(paste("Mean:", round(mean_val, 4)),
                  paste("Median:", round(median_val, 4)),
                  paste("SD:", round(sd_val, 4))),
       bty = "n",  # No box around legend
       cex = 1.1,
       text.col = "#333333")

dev.off()

end_time <- Sys.time()
execution_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

cat("Plot saved to:", output_file, "\n")
cat("Execution time:", execution_time, "seconds\n")

