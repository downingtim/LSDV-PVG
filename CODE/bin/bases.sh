#!/bin/bash

# Read the content of file_fa into a variable
read $content

# Remove 'N' characters
content_no_N=$(echo "$content" | tr -d 'N')

# Calculate the total length
total_length=$(echo "$content" | wc -c)
total_length=$((total_length - 10))

# Calculate the gap-free length
gap_free_length=$(echo "$content_no_N" | wc -c)
gap_free_length=$((gap_free_length - 10))

# Print the total length and gap-free length
echo "Total length  = $total_length" > bases.stats.txt
echo "Gap-free length  = $gap_free_length" >> bases.stats.txt

exit 1;
