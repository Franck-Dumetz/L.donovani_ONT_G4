# _L. donovani_ UTR mapping (using ONT Direct RNA sequencing)

BioProject SUB15103653

Table of content: <br />
[Basecalling ONT reads and mapping to the genome](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#basecalling-ont-reads-and-mapping-to-the-genome)<br />
[Isolating reads with L. donovani SL in 5'](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-reads-with-l-donovani-sl-in-5) <br />
  - [Isolating poly- and mono-cistronic reads](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-poly--and-mono-cistronic-reads) <br />
  - [Isolating reads with a SL](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-reads-with-a-sl) <br />
  - [Verify the SL is the same for all reads](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#verify-the-sl-is-the-same-for-all-reads) <br />

[Transcript evidence annotation and clean up](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#transcript-evidence-finding-and-clean-up) <br />
[Finding coding sequences](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#finding-coding-sequence) <br />
[Finding lncRNA in Leishmania genome](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#finding-lncRNA-in-leishmania-genome) <br />

Software requirements: <br />
- emboss-6.6.0 fuzznuc <br />
- gffread-0.12.7 <br />
- guppy-6.4.2 <br />
- feartureCounts-v2.0.6 <br />
- interproscan <br />
- manimap2.1 <br />
- meme-5.5.5 <br />
- ncbi-blast+-2.14.0 <br />
- samtools-1.20 <br />
- seqtk-1.0-r63-dirty <br />
- stringtie-2.2.1 <br />
- transdecoder-5.7.1 <br />

## Basecalling ONT reads and mapping to the genome

on a GPU <br />
[guppy.sh](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/guppy.sh) <br />
from the "pass" folder
```
cat *.fastq > experiment_name.fastq
```
map to the genome using minimap2 with [minimap2.sh](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/minimap.sh)<br />
```
samtools view -bhF 2308 sam_file.sam | samtools sort -o bam_file.bam
samtools index bam_file.bam
rm *.sam
```

## Isolating poly- and mono-cistronic reads

### Transfer of protein annoation from the reference genome into a new genome using gff

```
ref_prot=path_to_ref_genome/AnnotatedProteins.fasta
assembly=path_to_genome_assembly

/usr/local/packages/ncbi-blast+-2.14.0/bin/tblastn -query $ref_prot -subject $assembly -outfmt 6 -out ./tblastn.out
```
use [gff_from_blast.pl](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/gff_from_blast.pl) to transfer the tblastn result into a gff<br />


### Isolation of poly- and mono-cistronic

To isolate the polycistronic reads from the monocystronic ones, use the bam file previously generated and the gff generated above <br />
use [finding_polycistron.pl](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/finding_polycistron) <br />
it outputs 2 different sam files with the different reads and also the number for each category <br />
```
samtools view -bhF 2308 sam_file.sam | samtools sort -o bam_file.bam
samtools index bam_file.bam
rm *.sam
```
At this step you should have a polycistron.bam and a nomocistron.bam <br />

## Isolating reads with _L. donovani_ SL in 5'
Use the monocistron.bam <br /> 
Fuzznuc doesn't read bam files nor fastq files with U in it. <br />
Two different ways: <br />
### Work from the orginial fastq into fasta and change all Us into Ts
```
seqtk seq -a input.fastq > output.fasta

awk '/^[^>]/{ gsub(/U/,"T") }1' file.fasta >newfile.fasta
```
### Converting the bam file into a fastq file
```
bedtools bamtofastq -i monocistron.bam -fq monocistron.fastq
```
### Isolating reads with a SL
```
input=path/monocistron.fastq

fuzznuc -sequence $input -pattern ATAAGTATCAGTTTCTGTACTTTATTG -pmismatch 6 -outfile monocistron_SL.fuzznuc
```
The output needs to be refined to only have read names, start and end, start and end of motif and the motif found and searched for <br />
```
grep 'Sequence' --no-group-separator -A1 monocistron_SL.fuzznuc | grep -v 'HitCount' | grep -v 'Start' | awk '{printf "%s%s",$0,NR%2?"\t":RS}' > monocistron_SL.tsv
```
Then sort reads with a start of motif no futher than 6 nucleotides from the start <br />
Extract only read names from that file <br />
```
awk '{print $1}' monocistron_SL_filtered.txt |awk 'NR>1' > monocistron_SL_ReadName.txt
```
Create a new fastq file using the exported read names and the original monocistron.fastq <br />
```
grep --no-group-separator -A3 -f 'monocistron_SL_ReadName.txt' monocistron.fastq > monocistron_SL.fastq
```
Realign with [minimap2.sh](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/minimap.sh) <br />
and samtools (as previously described) <br />
creating monocistron_SL.bam <br />

### Verify the SL is the same for all reads
Isolating the 50 first nucleotides of each read <br />
```
seqkit subseq -r 1:50 monocistron_SL.fasta > 50first_nuc_SL.fasta
```
Using meme for motif search <br />
```
meme 50first_nuc_SL.fasta -dna -oc . -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun classic -markov_order 0
```
### Read lenght in the dataset
Longest read 
```
samtools view /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Fuzznuc/Ld1S_3ONT_monocistron.bam | awk '{print length($10)}' | sort -n | head -1
```
Shortest read 
```
samtools view /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Fuzznuc/Ld1S_3ONT_monocistron.bam | awk '{print length($10)}' | sort -nr | head -1
```
## Transcript evidence finding and clean up
```
input_bam=path/monocistron_SL.bam

stringtie-2.2.1/stringtie $input_bam --fr -f 0.5 -c 3 -l Ld1S -p 8 -L -R -o monocistron_SL.gtf
```
Data had to be manually clean up. Stringie is often confused with isoforms. Only kepp the longest transcript when it does not overlap with the neighboring reads. <br />
To create monocistron_SL_filtered.gtf <br />

## Finding coding sequence
```
gtf=path/monocistron_SL_filtered.gtf
assembly=path/assembly.fasta

transdecoder-5.7.1/util/gtf_genome_to_cdna_fasta.pl $gtf $assembly > Transdecoder_transcripts.fasta
```
```
transdecoder-5.7.1/TransDecoder.LongOrfs -t Transdecoder_transcripts.fasta
```
```
transdecoder-5.7.1/TransDecoder.Predict -t Transdecoder_transcripts.fasta --single_best_only
```
```
gtf=path/monocistron_SL_filtered.gtf

transdecoder-5.7.1/util/gtf_to_alignment_gff3.pl $gtf > monocistron_SL_filtered.gff3
```
```
transdecoder-5.7.1/util/cdna_alignment_orf_to_genome_orf.pl Transdecoder_transcripts.fasta.transdecoder.gff3 monocistron_SL_filtered.gff3 Transdecoder_transcripts.fasta > monocistron_SL.transdecoder.genome.gff3
```
TransDecoder outputs a lot of false positive. <br />
Transdcoder outputs 4 different files: .bed, .cds, .pep, .gff3 <br />
First, use the .pep file to compare with BPK282 using blatp: <br />
```
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastp -subject Transdecoder_transcripts.fasta.transdecoder.pep -query TriTrypDB-66_LdonovaniBPK282A1_AnnotatedProteins.fasta -outfmt 6 -out blastp.out
```
Then only keep all lines with a match value >95% in a blastp_filtered.out <br />
Isolate the headers from the .pep file <br />
```
awk '/^>/ { print $0 > "pep-headers.txt"; next } { print $0 > "pep-sequences.txt" }' Transdecoder_transcripts.fasta.transdecoder.pep

awk '{print S1}' blast_filtered.out > blast_LdName.txt
```
Use [Match_pep2blastp.py](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/match_pep2blastp.py) <br /> 
It outputs a file called complete_1Sfrom-282.txt that contains all BPK282 transfered to Ld1S (6201 protein transfered) <br />

Asign header to sequence using [pep_file_Ld1S.py](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/pep_file_Ld1S.py) <br />
Output file Ld1S_pep_from282.fasta contains all pep sequence transfered from BPK282 <br />
<br />
<br />
<br />
<br />
Isolating only transcripts annotation from the Stringtie output. This will remove all mention of UTRs <br />
```
awk '$3 == "gene" || $3 == "mRNA" || $3 == "exon" || $3 == "CDS" || $3 == ""' Ld1S_stg_filtered.transdecoder.genome.gff3 > Ld1S_stg_filtered.transdecoder.genome_short.gff3

```
Making a new gff3 with all BPK282 transdecoder detected annotations <br />
to only select BPK282 transferred annotation use [new_gff3_TransDecoderVS282.py](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/new_gff3_TransDecoderVS282.py) <br />
Clean up to remove all non CDS annotation with: <br />
```
grep -v 3prime Cleaned_up_transdecoderVSBPK282.gff3 >  Cleaned_up_transdecoderVSBPK282_1.gff3
grep -v 5prime Cleaned_up_transdecoderVSBPK282.gff3 >  Cleaned_up_transdecoderVSBPK282_2.gff3
grep -v -f blast_LdName.txt Cleaned_up_transdecoderVSBPK282_2.gff3 | head
```

## UTR position and length

Use the original Stringtie output with only the "transcript" lines and the gff3 files from TransDecoder where all the lines with "ORF type complete" are isolated. <br />
Then use [UTR_position-length](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/UTR_position-length.pl) <br />



## Filtering out reads with less than 5x coverage

```
samtools depth Ld1S_Ama_monocistron_SL.bam > Ld1S_Ama_monocistron_SL_coverage.txt
awk '$3 >=5' Ld1S_Ama_monocistron_SL_coverage.txt > Ld1S_Ama_monocistron_SL_FiltCoverage.txt
awk '{print $1"\t"$2-1"\t"$2}' Ld1S_Ama_monocistron_SL_FiltCoverage.txt > Ld1S_Ama_monocistron_SL_FiltCoverage.bed
bedtools intersect -abam Ld1S_Ama_monocistron_SL.bam -b Ld1S_Ama_monocistron_SL_FiltCoverage.bed > Ld1S_Ama_monocistron_SL_FiltCoverage.bam
samtools index Ld1S_Ama_monocistron_SL_FiltCoverage.bam
```

## Naming harmonisazion between Stringtie promastigotes and Stringtie amastigotes

To match names between the different files, we used [gtf_name_transfer.pl](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/gtf_name_transfer.pl) <br />
It will transfer the name of the promastigote stringtie onto the amastigote stringie <br />


## Finding lncRNA in Leishmania genome

## Identifying nucleic acid motif at splicing junction

Extracting 5'UTR bounderies to output chromosome name and ending position of the 5' UTR
```
grep 'exon' Ld1S_annotation_filtered_final_final.gtf | awk '{print $1"\t"$4}' > 5bounderies.txt
```
Extracting 3'UTR bounderies to output chromosome name and ending position of the 3' UTR
```
grep 'exon' Ld1S_annotation_filtered_final_final.gtf | awk '{print $1"\t"$5}' > 3bounderies.txt
```
Extracting 25 nucleotides before and after the splicing junction from the fasta genome
```
awk '{print $1 ":" $2-25 "-" $2+25}' 5bounderies.txt | xargs -I {} samtools faidx /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_genome/Flye_scaffold/Ld1S_assembly_final.fasta {} > 5bounderies_seq.fasta
```
Meme to find motifs 
```
/usr/local/packages/meme-5.5.5/bin/meme /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_UTR/5bounderies_seq.fasta -dna -oc . -mod zoops -nmotifs 10 -maxw 25 -objfun classic -markov_order 0
```
## Gene expression
```
featureCounts -s 2 -T 8 -a /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step4_final_removal/Ld1S_annotation_filtered_final_final.gtf -o Ld_3ONT_count.txt /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Fuzznuc/Ld1S_3ONT_monocistron.bam
```
## Gene ontology analysis
Prepare the data
```
grep -A1 --no-group-separator -E '\.1\.p1([^0-9]|$)' linear_Ld1S_pep.fasta |sed 's/\*$//' > onlyp1.txt
```
Use interproscan to find GO terms associated to each protein
```
interproscan.sh -i AMA_pep_clean.fasta -f TSV -goterms -pa -cpu 16 -o AMA_blast_resultsGO.tsv
```
```
cut -f1,14 AMA_blast_resultsGO2.tsv | grep "GO:" | less
```
Use topGO in R like in [topGo]()
