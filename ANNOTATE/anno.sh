#!/bin/bash

# Genome Annotation Processing Script
# Purpose: Process prokka-generated GFF files and create graph node annotations
# Author: [Your Name]
# Date: $(date +%Y-%m-%d)

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
REFERENCE_ID="KX894508"
INPUT_GFF="${REFERENCE_ID}.gff"
GRAPH_FILE="LSDV3.odgi"

# Output files
CLEAN_GFF="${REFERENCE_ID}_clean.gff"
ANNOTATIONS_BED="annotations.bed"
UNTANGLE_BED="${REFERENCE_ID}_untangle.bed"
INTERSECT_BED="nodes_with_annotations.bed"
FINAL_OUTPUT="graph_node_annotations.tsv"

# ============================================================================
# FUNCTIONS
# ============================================================================

log_step() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_file_exists() {
    if [[ ! -f "$1" ]]; then
        echo "Error: Required file '$1' not found" >&2
        exit 1
    fi
}

# ============================================================================
# MAIN PROCESSING
# ============================================================================

log_step "Starting genome annotation processing pipeline"

# Load required modules
log_step "Loading required modules"
module load tools
module load prokka

# Check input files exist
log_step "Checking input files"
check_file_exists "$INPUT_GFF"
check_file_exists "$GRAPH_FILE"

# Step 1: Clean GFF file
log_step "Cleaning GFF file: $INPUT_GFF -> $CLEAN_GFF"
perl -pe 's/gnl\|Prokka\|KCBEBCIG_1/KX894508/g' "$INPUT_GFF" | \
    sed 's/^>//' > "$CLEAN_GFF"

# Step 2: Convert GFF to BED format
log_step "Converting GFF to BED format: $CLEAN_GFF -> $ANNOTATIONS_BED"
awk '
BEGIN { OFS="\t" }
$1 !~ /^#/ {
    start = ($4 > 0) ? $4 - 1 : 0
    end = $5
    if (start >= 0 && end > start) {
        print $1, start, end, $9
    }
}' "$CLEAN_GFF" > "$ANNOTATIONS_BED"

# Step 3: Extract untangled nodes from graph
log_step "Extracting untangled nodes: $GRAPH_FILE -> $UNTANGLE_BED"
odgi untangle -i "$GRAPH_FILE" -r "$REFERENCE_ID" > "$UNTANGLE_BED"

# Step 4: Intersect annotations with untangled nodes
log_step "Intersecting annotations with nodes: $UNTANGLE_BED + $ANNOTATIONS_BED -> $INTERSECT_BED"
bedtools intersect \
    -a "$UNTANGLE_BED" \
    -b "$ANNOTATIONS_BED" \
    -wa -wb > "$INTERSECT_BED"

# Step 5: Parse annotations and create final output
log_step "Parsing annotations and creating final output: $INTERSECT_BED -> $FINAL_OUTPUT"
awk '
BEGIN { 
    OFS="\t"
    print "node_id", "id", "locus_tag", "gene", "product"
}
{
    node_id = $4
    attr = $14
    
    # Initialize variables
    id = ""
    locus = ""
    gene = ""
    product = ""
    
    # Split attributes by semicolon
    n = split(attr, arr, ";")
    
    # Parse each attribute
    for (i = 1; i <= n; i++) {
        if (arr[i] ~ /^ID=/) {
            id = arr[i]
            sub(/^ID=/, "", id)
        }
        else if (arr[i] ~ /^locus_tag=/) {
            locus = arr[i]
            sub(/^locus_tag=/, "", locus)
        }
        else if (arr[i] ~ /^gene=/) {
            gene = arr[i]
            sub(/^gene=/, "", gene)
        }
        else if (arr[i] ~ /^product=/) {
            product = arr[i]
            sub(/^product=/, "", product)
        }
    }
    
    # Only output gene entries
    if (id ~ /_gene$/) {
        print node_id, id, locus, gene, product
    }
}' "$INTERSECT_BED" > "$FINAL_OUTPUT"

echo "Summary of outputs:"
echo "  - Clean GFF file: $CLEAN_GFF"
echo "  - BED annotations: $ANNOTATIONS_BED"
echo "  - Untangled nodes: $UNTANGLE_BED"
echo "  - Intersected data: $INTERSECT_BED"
echo "  - Final annotations: $FINAL_OUTPUT"
