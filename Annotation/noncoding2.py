import re

# ---- STEP 1 ----
# Load ncRNA transcript IDs

ncRNA_file = "/Users/fdumetz/Desktop/Desktop - fdumetzm4-osx/leish_gff/ncRNA_transcript_ids.txt"
ncRNA_ids = set()

with open(ncRNA_file, "r") as f:
    for line in f:
        tid = line.strip()
        if tid:
            ncRNA_ids.add(tid)

print(f"Loaded {len(ncRNA_ids)} transcript IDs from ncRNA_transcript_ids.txt")

# ---- STEP 2 ----
# Scan GFF and extract exons whose Parent is in our ncRNA list

gff_file = "/Users/fdumetz/Desktop/Desktop - fdumetzm4-osx/leish_gff/Ld1S_full_ortho_70_90.gff"
output_file = "/Users/fdumetz/Desktop/Desktop - fdumetzm4-osx/leish_gff/exon_parents_for_ncRNA.txt"

count = 0
exon_lines_checked = 0

with open(gff_file, "r") as f, open(output_file, "w") as out:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue

        feature = fields[2]
        if feature != "exon":
            continue

        exon_lines_checked += 1

        attr_field = fields[8]

        # Find Parent attribute
        match = re.search(r'Parent=([^;]+)', attr_field)
        if match:
            parent_id = match.group(1).strip()

            # Remove prefixes if present (e.g. transcript:)
            parent_no_prefix = re.sub(r'^\w+:', '', parent_id)

            # Remove trailing .p# if present
            parent_clean = re.sub(r'\.p\d+$', '', parent_no_prefix)

            # Debug for first few lines
            if exon_lines_checked <= 10:
                print(f"Line #{exon_lines_checked}")
                print(f"Original Parent: {parent_id}")
                print(f"Cleaned Parent: {parent_clean}")

            if parent_clean in ncRNA_ids:
                out.write(parent_id + "\n")
                count += 1
                print(f"*** MATCHED: original Parent saved: {parent_id}")

print(f"Total exon lines scanned: {exon_lines_checked}")
print(f"Total matching Parent IDs written to output: {count}")
print(f"Output saved to: {output_file}")