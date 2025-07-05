#!/usr/bin/env python3

import sys

if len(sys.argv) != 4:
    print("Usage: python add_attribute2gff.py attibutes_matchingg.tsv input.gff output.gff")
    sys.exit(1)

orthologs_file = sys.argv[1]
gff_file = sys.argv[2]
output_file = sys.argv[3]

# Load ortholog mapping
orthologs = {}
with open(orthologs_file) as f:
    for line in f:
        if line.strip() == "":
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            leish_id, bpk_id = parts[0], parts[1]
            orthologs[leish_id] = bpk_id

with open(gff_file) as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue

        cols = line.strip().split("\t")
        if len(cols) < 9:
            outfile.write(line)
            continue

        attributes = cols[8]

        # Check if any known gene ID matches this line
        found = False
        for leish_id in orthologs:
            if f"ID={leish_id}" in attributes:
                attributes += f";LdBPK_ortho={orthologs[leish_id]}". ### change for any attibute name 
                found = True
                break

        cols[8] = attributes
        outfile.write("\t".join(cols) + "\n")
