#!/usr/bin/env python3

# File paths for the input files and output file
gff3_file_path = '/PATH/Ld1S_stg_filtered.transdecoder.genome_short.gff3'
blast_file_path = '/PATH/blast_LdName.txt'
output_gff3_file = 'extracted_gff3_records.gff3'

# Extracting the IDs from blast_LdName.txt
blast_ids = set()

with open(blast_file_path, 'r') as blast_file:
    for line in blast_file:
        blast_ids.add(line.strip())  # Add each line as an ID without trailing characters

# Reading the GFF3 file and writing lines that contain any of the blast IDs
with open(gff3_file_path, 'r') as gff3_file, open(output_gff3_file, 'w') as output_file:
    for line in gff3_file:
        # Check if any of the blast IDs are in the line
        if any(blast_id in line for blast_id in blast_ids):
            output_file.write(line)

print(f"Extracted GFF3 records saved to {output_gff3_file}")

