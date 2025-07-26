#!/usr/bin/env python3

import sys

def calculate_genome_lengths(input_file, output_file="genome.lengths.txt"):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        sequences = infile.read().split(">")
        for seq in sequences[1:]:  # Skip the first empty split
            header, *sequence = seq.split("\n", 1)
            sequence = "".join(sequence).replace("\n", "")
            outfile.write(f"{header}\t{len(sequence)}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_fasta>")
    else:
        calculate_genome_lengths(sys.argv[1])
