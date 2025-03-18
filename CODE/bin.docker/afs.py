#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

vcf_file = "gfavariants.vcf"
total_samples = int(sys.argv[2])
output_file = "out.vcf.txt"

# Initialize dictionary to store counts
h = {}

# Read input VCF file
with open(vcf_file, "r") as infile:
    lines = infile.readlines()

# Process lines (skipping last one)
for line in lines[:-1]:  
    if "snv" in line and not line.startswith("#"):  # Ignore comment lines
        fields = line.split()
        pos = int(fields[1])  # POS field
        
        if pos > 300:  # Apply filter
            h[pos] = h.get(pos, 0) + 1

# Write to output file and store frequencies for plotting
frequencies = []
with open(output_file, "w") as outfile:
    for pos in sorted(h.keys()):
        freq = h[pos] / total_samples
        frequencies.append(freq)
        outfile.write(f"{freq}\n")

print(f"Output written to {output_file}")

# Plot the frequency distribution
plt.figure(figsize=(6, 4))
plt.hist(frequencies, bins=100, range=(0,1), edgecolor='black', alpha=0.75)
plt.xlabel("Frequency")
plt.ylabel("Number of observations")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.savefig("frequency_distribution.png", dpi=400)
plt.show()
