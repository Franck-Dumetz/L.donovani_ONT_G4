#!/usr/bin/env python3

# File paths for the input files
pep_file_path = '/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/linear_Ld1S_Transdecoder_transcripts.fasta.transdecoder.pep'
complete_file_path = '/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/complete_1Sfrom282.txt'
output_sequences_file = 'Ld1S_pep_from282.fasta'

# Extracting the IDs from the complete_1Sfrom-282.txt
complete_ids = set()

with open(complete_file_path, 'r') as complete_file:
    for line in complete_file:
        if line.startswith('>'):
            complete_ids.add(line.split()[0][1:])  # Extract ID without '>' character

# Reading through the pep file and extracting sequences with matching IDs
with open(pep_file_path, 'r') as pep_file, open(output_sequences_file, 'w') as output_file:
    write_sequence = False
    for line in pep_file:
        if line.startswith('>'):
            sequence_id = line.split()[0][1:]  # Extract ID without '>' character
            if sequence_id in complete_ids:
                write_sequence = True
                output_file.write(line)
            else:
                write_sequence = False
        elif write_sequence:
            output_file.write(line)

print(f"Extracted sequences saved to {output_sequences_file}")

