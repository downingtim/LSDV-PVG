import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_file> <divisor>")
    sys.exit(1)

# Get input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[1] + ".out"

# Get the divisor from command-line arguments
divisor = int(sys.argv[2])

# Open input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Initialize a dictionary to store counts
    h = {}

    # Read all lines from the input file
    lines = infile.readlines()

    # Loop through each line
    for line in lines:
        # Split the line into fields
        fields = line.split()

        # Check if the line contains 'snv'
        if 'snv' in fields:
            # Extract relevant fields
            chrom, pos = fields[0], int(fields[1])

            # Check if position is within range
            if 1000 < pos < 148001:
                # Increment count for the position
                h[pos] = h.get(pos, 0) + 1

    # Sort dictionary keys and write counts to output file
    for key in sorted(h.keys()):
        if h[key] > 0:
            outfile.write(f"{h[key] / divisor}\n")

sys.exit(0)
