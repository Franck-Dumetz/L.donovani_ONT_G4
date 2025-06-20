#!/usr/bin/env python3

### Matching pep headers to blastp output

# File paths for the input files
blast_file_path = '/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/blast_LdName.txt'
pep_headers_path = '/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/pep-headers.txt'
output_file_path = '/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/complete_1Sfrom-282.txt'  # Output file


# Load the blast IDs from blast_LdName.txt into a set
barcode_dict = {}

blast_file = open(blast_file_path, 'r')
for line in blast_file:
    barcode_dict[line.strip()] = 1

# Filter the pep-headers based on blast names and presence of the word 'complete', and save it to a new file
with open(pep_headers_path, 'r') as pep_file, open(output_file_path, 'w') as output_file:
    for line in pep_file:
        split_line = line.split(" ")
        barcode = split_line[0][1:]

        # Check if any of the blast IDs are in the line and if 'complete' is in the line
        if barcode in barcode_dict and 'complete' in line:
            output_file.write(line)

print(f"Filtered results saved to {output_file_path}")

