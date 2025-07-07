#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict
import re
import os
import sys

def parse_gff3(gff3_file, feature_type):
    transcript_exons = defaultdict(list)
    with open(gff3_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2].lower() != feature_type.lower():
                continue
            chrom, _, _, start, end, _, strand, _, attr = fields
            start, end = int(start), int(end)

            # Try transcript_id, Parent, or ID
            tid_match = re.search(r'transcript_id[ ="]*([^";\s]+)', attr)
            if not tid_match:
                tid_match = re.search(r'Parent=([^";\s]+)', attr)
            if not tid_match:
                tid_match = re.search(r'ID=([^";\s]+)', attr)
            if tid_match:
                tid = tid_match.group(1)
                transcript_exons[tid].append((chrom, start, end))
    return transcript_exons

def process_bams(transcript_exons, bam_files, feature_type, allow_multimapping):
    all_counts = defaultdict(lambda: defaultdict(int))
    sample_names = [os.path.basename(bam).replace(".bam","") for bam in bam_files]

    for sample_name, bam_path in zip(sample_names, bam_files):
        print(f"Processing {sample_name}...")

        bam = pysam.AlignmentFile(bam_path, "rb")

        if allow_multimapping:
            # ORIGINAL BEHAVIOR
            for tid, exons in transcript_exons.items():
                seen_reads = set()
                for chrom, start, end in exons:
                    try:
                        for read in bam.fetch(chrom, start - 1, end):
                            if read.is_unmapped:
                                continue
                            if read.query_name in seen_reads:
                                continue
                            if read.reference_start < end and read.reference_end > start:
                                all_counts[tid][sample_name] += 1
                                seen_reads.add(read.query_name)
                    except ValueError:
                        continue
        else:
            # BLOCK MULTIMAPPING ACROSS FEATURES
            read_to_tids = defaultdict(set)

            for tid, exons in transcript_exons.items():
                for chrom, start, end in exons:
                    try:
                        for read in bam.fetch(chrom, start - 1, end):
                            if read.is_unmapped:
                                continue
                            if read.reference_start < end and read.reference_end > start:
                                read_to_tids[read.query_name].add(tid)
                    except ValueError:
                        continue

            multimapped_reads = 0
            counted_reads = set()

            for read_name, tid_set in read_to_tids.items():
                if len(tid_set) == 1:
                    tid = next(iter(tid_set))
                    if read_name not in counted_reads:
                        all_counts[tid][sample_name] += 1
                        counted_reads.add(read_name)
                else:
                    multimapped_reads += 1

            print(f"[{sample_name}] Reads discarded due to overlapping multiple {feature_type} features: {multimapped_reads}")

        bam.close()

    return all_counts, sample_names

def write_output(all_counts, transcript_exons, sample_names, output_file):
    # compute total reads per sample
    total_counts = {sample: 0 for sample in sample_names}
    for tid_counts in all_counts.values():
        for sample, count in tid_counts.items():
            total_counts[sample] += count

    # write file
    with open(output_file, "w") as out:
        out.write("transcript_id\t" +
                  "\t".join(f"{s}_raw" for s in sample_names) + "\t" +
                  "\t".join(f"{s}_CPM" for s in sample_names) + "\n")

        for tid in sorted(transcript_exons):
            raw_counts = [all_counts[tid][s] for s in sample_names]
            cpm_values = [
                f"{(raw / total_counts[sample]) * 1e6:.2f}" if total_counts[sample] > 0 else "0.00"
                for raw, sample in zip(raw_counts, sample_names)
            ]
            raw_str = "\t".join(str(rc) for rc in raw_counts)
            cpm_str = "\t".join(cpm_values)
            out.write(f"{tid}\t{raw_str}\t{cpm_str}\n")

    print(f"Output written to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Count reads overlapping genomic features from a GFF3 file across multiple BAMs."
    )
    parser.add_argument("-f", "--feature", required=True, help="Feature type to extract (e.g. exon, CDS, gene)")
    parser.add_argument("-b", "--bam", nargs='+', required=True, help="One or more BAM files")
    parser.add_argument("-g", "--gff3", required=True, help="GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("--allow-multimapping", action="store_true",
                        help="Allow reads to be counted in multiple features. If absent, multimapping reads are discarded.")

    args = parser.parse_args()

    if not os.path.exists(args.gff3):
        print(f"ERROR: GFF3 file not found: {args.gff3}")
        sys.exit(1)

    for bam_path in args.bam:
        if not os.path.exists(bam_path):
            print(f"ERROR: BAM file not found: {bam_path}")
            sys.exit(1)

    transcript_exons = parse_gff3(args.gff3, args.feature)
    print(f"Found {len(transcript_exons)} {args.feature} features grouped by transcript.")

    all_counts, sample_names = process_bams(transcript_exons, args.bam, args.feature, args.allow_multimapping)
    write_output(all_counts, transcript_exons, sample_names, args.output)

if __name__ == "__main__":
    main()
