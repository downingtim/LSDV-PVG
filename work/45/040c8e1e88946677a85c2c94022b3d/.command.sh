#!/bin/bash -ue
bgzip -@ 4 communities.genomes.fasta
    samtools faidx communities.genomes.fasta.gz # index

    # Community detection based on p/w alignments based on 90% ID at mash level with 6
    # mappings per segment, using k-mer of 19 and window of 67
    wfmash communities.genomes.fasta.gz -p 90 -n 6 -t 44 -m > genomes.mapping.paf

    # Convert PAF mappings into a network:
    paf2net.py -p genomes.mapping.paf

    # make plot
    net2communities.py -e genomes.mapping.paf.edges.list.txt -w genomes.mapping.paf.edges.weights.txt -n genomes.mapping.paf.vertices.id2name.txt --plot

    # mash-based partitioning
    mash dist communities.genomes.fasta.gz communities.genomes.fasta.gz -s 10000 -i > genomes.distances.tsv

    # get distances
    mash2net.py -m genomes.distances.tsv

    output_file="genomes.communities.mash.paths.txt"

    # Initialize a flag to check if any community files were found
    found_community_files=0

    for i in {0..3}; do
    # Check if the community file exists
        if [[ -f "genomes.mapping.paf.edges.weights.txt.community.${i}.txt" ]]; then
    # Extract unique chromosome names from the community file
            chromosomes=$(cut -f 3 -d '#' "genomes.mapping.paf.edges.weights.txt.community.${i}.txt" | sort | uniq | tr '
' ' ')
            echo "community ${i} --> ${chromosomes}"
            found_community_files=1
        else
            echo "community ${i} --> Community file not found"
        fi
    done

    # Check if any community files were found
    if [[ $found_community_files -eq 0 ]]; then
        echo "No community files found. Exiting."
    fi  
    # Redirect the output to the specified file
    echo "Output written to $output_file"
    # community_numbers.sh end

    pafgnostic --paf genomes.mapping.paf > genomes.mapping.paf.txt 
    #  pafgnostic outputs a tab-separated table directly to the console with the following columns:
    # - Query name, start, end, strand
    # - Target name, start, end
    # - Number of matches, mismatches, insertions, deletions
    # - Unique matches, mismatches, insertions, deletions
    # - Longest, shortest, and average lengths of matches, mismatches, insertions, deletions
