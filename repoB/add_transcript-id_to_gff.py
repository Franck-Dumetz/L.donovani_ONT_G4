input_file = '/Users/fdumetz/Desktop/Leish_lncRNA/Ld1S_stg_filtered.transdecoder.genome.gff3'
output_file = '/Users/fdumetz/Desktop/Leish_lncRNA/Ld1S_stg_filtered.transdecoder.genome.fixed.gff3'

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith("#"):
            # Write comment lines as is
            outfile.write(line)
        elif "ID=" in line:
            parts = line.strip().split("\t")
            attributes = parts[-1]  # Extract attributes column
            if "transcript_id=" not in attributes:
                # Extract the ID and add as transcript_id
                id_attr = [attr for attr in attributes.split(";") if attr.startswith("ID=")]
                if id_attr:
                    transcript_id = id_attr[0].replace("ID=", "")
                    attributes += f";transcript_id={transcript_id}"
                    parts[-1] = attributes
            outfile.write("\t".join(parts) + "\n")
        else:
            # Write other lines as is
            outfile.write(line)

print(f"Fixed GFF3 written to {output_file}")
