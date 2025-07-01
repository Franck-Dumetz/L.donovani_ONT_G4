#!/usr/bin/env python3

# trnascan_to_gff.py
# Convert tRNAscan-SE tabular output to GFF3

import sys

# Usage: python trnascan_to_gff.py trnascan_output.txt output.gff3

if len(sys.argv) != 3:
    print("Usage: python trnascan_to_gff.py <trnascan_output.txt> <output.gff3>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as infile, open(output_file, "w") as out:
    out.write("##gff-version 3\n")
    
    for line in infile:
        # Skip headers or empty lines
        if line.strip() == "" or line.startswith("-") or line.startswith("Sequence"):
            continue

        cols = line.strip().split()
        
        seqid = cols[0]
        trna_num = cols[1]
        start = cols[2]
        end = cols[3]
        trna_type = cols[4]
        codon = cols[5]
        anticodon = cols[6]
        intron = cols[7] if len(cols) > 7 else "0"

        strand = "+"
        if int(start) > int(end):
            start, end = end, start
            strand = "-"
        
        attributes = f"ID=tRNA{trna_num};Product=tRNA-{trna_type};Anticodon={anticodon};Codon={codon};Intron={intron}"
        
        gff_line = f"{seqid}\ttRNAscan-SE\ttRNA\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
        
        out.write(gff_line)

print(f"Output written to {output_file}")