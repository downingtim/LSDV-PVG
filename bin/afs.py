#!/usr/bin/env python3

import sys

def calculate_afs(input_file, genome_length, output_file=None):
    if output_file is None:
        output_file = f"{input_file}.out"
    
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        data = infile.readlines()
        afs_counts = {}
        genome_end_threshold = genome_length - 100
        
        for line in data:
            if "snv" in line:
                fields = line.split()
                pos = int(fields[1])
                
                if 1000 < pos < genome_end_threshold:
                    afs_counts[pos] = afs_counts.get(pos, 0) + 1
        
        for pos in sorted(afs_counts.keys()):
            outfile.write(f"{afs_counts[pos] / int(genome_length)}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_fasta> <genome_length>")
    else:
        calculate_afs(sys.argv[1], int(sys.argv[2]))
