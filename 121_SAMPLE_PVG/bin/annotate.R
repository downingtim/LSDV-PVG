#!/usr/bin/Rscript

library(rentrez) # Load the rentrez package

# Function to fetch GenBank file for a given virus
fetch_genbank_file <- function(virus) {
  # Define the output file name
  output_file <- paste0("../../../CURRENT/PANGROWTH/SEQS/",virus, ".gb")

  # Define the accession number based on the virus
  accession_number <- switch(virus,
                             "lsdv" = "KX894508",
                             "gpv"  = "NC_004003.1",
                             "sppv" = "NC_004002.1",
                             stop("Invalid virus argument"))

  # Specify the database (nuccore for nucleotide sequences) and rettype (GenBank format)
  database <- "nuccore"
  rettype <- "gb"
  # Construct the query to fetch the GenBank file
  query <- paste(accession_number, "txid10241[ORGN]", sep = " AND ")
  # Fetch the GenBank file
  genbank_data <- entrez_fetch(db = database, id = query, rettype = rettype, mode = "text")
  # Write the GenBank data to a file
  writeLines(genbank_data, con = output_file)
  cat("GenBank file", output_file, "downloaded successfully.\n")
}

# Function to run Prokka for each .fa file
run_prokka <- function(file) {
  c <- basename(tools::file_path_sans_ext(file))
  # Append "_PROKKA" to the modified file name
  b <- paste0("../../../CURRENT/PANGROWTH/SEQS/",c, "_PROKKA")
  
  # Print the original file name
  cat(file, "\n")
  
  # Run prokka command
  print(paste('prokka ',' --kingdom Viruses --proteins KX894508.gb --gffver 3 --usegenus --outdir ',b,' --force --genus Capripoxvirus --prefix ',c,' ',file, sep=""))
  system2(command="prokka", args=c(' --kingdom Viruses --proteins KX894508.gb --gffver 3 --usegenus --outdir ',b,' --force --genus Capripoxvirus --prefix ',c,' ',file), wait = T, stdout = paste0(b, ".out"), stderr = paste0(b, ".err"))
}

# Check if the correct number of arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 1) {
  stop("Usage: ./combined_script.R [lsdv | gpv | sppv]") }

# Get the virus argument from the command line
virus_arg <- tolower(commandArgs(trailingOnly = TRUE)[1])

# Check if the virus argument is LSDV, GPV, or SPPV and fetch the corresponding GenBank file
fetch_genbank_file(virus_arg)

# Set the working directory to the parent directory of "bin"
fa_files <- list.files("../../../CURRENT/PANGROWTH/SEQS", pattern = "\\.fa$", full.names = TRUE)

# Iterate over each .fa file and run Prokka
for (file in fa_files) {  run_prokka(file) }

q()
y