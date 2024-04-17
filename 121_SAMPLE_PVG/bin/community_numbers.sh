#!/bin/bash

# Define the output file path
output_file="../../../CURRENT/COMMUNITIES/genomes.communities.mash.paths.txt"

# Initialize a flag to check if any community files were found
found_community_files=0

# Iterate over the range of community file indices, e.g., 0 to 3
for i in {0..3}; do
    # Check if the community file exists
    if [[ -f "../../../CURRENT/COMMUNITIES/genomes.mapping.paf.edges.weights.txt.community.${i}.txt" ]]; then
        # Extract unique chromosome names from the community file
        chromosomes=$(cut -f 3 -d '#' "../../../CURRENT/COMMUNITIES/genomes.mapping.paf.edges.weights.txt.community.${i}.txt" | sort | uniq | tr '\n' ' ')
        echo "community $i --> $chromosomes"
        found_community_files=1
    else
        echo "community $i --> Community file not found"
    fi
done

# Check if any community files were found
if [[ $found_community_files -eq 0 ]]; then
    echo "No community files found. Exiting."
    exit 1
fi

# Redirect the output to the specified file
echo "Output written to $output_file"
exit 0
