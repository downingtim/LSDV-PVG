#!/usr/bin/env python3

import sys

def calculate_base_lengths(input_file, output_file="../../../CURRENT/bases.txt"):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        data = infile.read()
        total_length = len(data) - 10
        gap_free_length = len(data.replace("N", "")) - 10
        outfile.write(f"Total length  = {total_length}\tGap-free length  = {gap_free_length}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_fasta>")
    else:
        calculate_base_lengths(sys.argv[1])
