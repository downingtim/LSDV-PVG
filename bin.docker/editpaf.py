#!/usr/bin/env python3
import sys
input_paf = sys.argv[1]
bed_file = sys.argv[2]

Small_Genomes = {}

with open(input_paf, 'r') as infile, open(bed_file, 'r') as bedfile:
        for line in bedfile:
            # Skip header lines
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            Small_Genomes[fields[0]] = int(fields[2])  # Store the length of each genome

        for line in infile:
            # Skip header lines
            if line.startswith('#'):
                outfile.write(line)
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 12:  # PAF requires at least 12 fields
                continue
                
            # Parse PAF fields
            query_name = fields[0]
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_length = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            sketch_match = int(fields[9])
            alignment_block_from_query = int(fields[10])
            mapping_quality = int(fields[11])

            if target_name in Small_Genomes:
                target_actual_length = Small_Genomes[target_name]
                if target_start > target_actual_length: continue
                if target_end > target_actual_length: target_end = target_actual_length
                target_length = target_actual_length
                    
            if query_name in Small_Genomes:
                query_actual_length = Small_Genomes[query_name]
                if query_start > query_actual_length: continue
                if query_end > query_actual_length: query_end = query_actual_length
                query_length = query_actual_length
                if alignment_block_from_query > query_actual_length:
                    alignment_block_from_query = query_actual_length
            
            print('\t'.join([
                query_name, str(query_length), str(query_start), str(query_end), strand,
                target_name, str(target_length), str(target_start), str(target_end),
                str(sketch_match), str(alignment_block_from_query), str(mapping_quality),
                fields[12], fields[13]
            ]))
