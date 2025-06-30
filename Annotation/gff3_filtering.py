keep_ids = set(
    line.strip()
    for line in open("/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/250630_Ld1S_names_from282.txt")
)

with open("/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/Ld1S_stg_filtered.transdecoder.genome.gff3") as infile, \
     open("/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/Ld1S_stg_filtered_transdecoder.gff3", "w") as outfile:
    
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue

        # Split attributes field (9th column) to look for IDs
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]

        found = False
        for keep_id in keep_ids:
            if keep_id in attributes:
                found = True
                break

        if found:
            outfile.write(line)
