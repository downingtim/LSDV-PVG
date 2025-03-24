#!/usr/bin/env python
import re
import statistics
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def extract_year(feature, qualifiers):
    if feature in qualifiers and qualifiers[feature]:
        match = re.search(r'\b(19\d{2}|20\d{2})\b', qualifiers[feature][0])
        return match.group(0) if match else None
    return None

# Initialize missing lists
RC_List = []  # Define this properly if needed
Lengths = []

# Parse the GenBank file
with open(sys.argv[1], "r") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        Collection_Date = ""
        Host = ""

        for feature in record.features:
            if 'host' in feature.qualifiers:
                Host = feature.qualifiers['host'][0]
            if 'collection_date' in feature.qualifiers:
                Collection_Date = feature.qualifiers['collection_date'][0]
            else:
                for key in ["collected_by", "note", "isolate"]:
                    extracted_date = extract_year(key, feature.qualifiers)
                    if extracted_date:
                        Collection_Date = extracted_date
                        break

        if record.id in RC_List:
            record.seq = record.seq.reverse_complement()

        seq_length = len(record.seq)
        if 4000 <= seq_length <= 1470000:
            Lengths.append(seq_length)
            print(f">{record.id} {record.description}_{Collection_Date}_{Host}")

            # Ensure sequence is long enough before slicing
            if seq_length > 2002:
                print(record.seq.lower()[50:-50] + 'n' * 50)
            else:
                print(record.seq.lower())  # Print full sequence if too short
