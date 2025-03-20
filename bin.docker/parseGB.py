#!/usr/bin/env python
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq


def extract_year(feature,features):
    if feature in features:
       match = re.search(r'\b(19\d{2}|20\d{2})\b',  features[feature][0])
       return match.group(0) if match else None
    else:
       return(None)

RC_List = ["OK422492.1","OK422493.1","ON400507.1","PP145891.1","OM373209.1"]

# Parse the GenBank file


with open(sys.argv[1], "r") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        Collection_Date = ""; Host="" 
        for feature in record.features:
            if 'host' in feature.qualifiers:
                Host = feature.qualifiers['host'][0]
            if 'collection_date' in feature.qualifiers:
                Collection_Date = feature.qualifiers['collection_date'][0]
            else:
                for Key in ["collected_by","note","isolate"]:
                    C = extract_year(Key,feature.qualifiers)
                    if C != None: 
                        Collection_Date = C
                        break
                    
        if record.id in RC_List: record.seq = record.seq.reverse_complement()
        if len(record.seq) <= 180000 and len(record.seq) >= 147000:
            print(">%s %s_%s"%(record.id,record.description,"_".join([Collection_Date,Host])))
            print(format('n'*1000 + record.seq.lower()[1000:-1001] + 'n'*1001))






