# _L. donovani_ UTR mapping (using ONT Direct RNA sequencing)

Table of content: <br />
[Basecalling ONT reads and mapping to the genome](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#basecalling-ont-reads-and-mapping-to-the-genome)<br />
[Isolating reads with L. donovani SL in 5'](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-reads-with-l-donovani-sl-in-5) <br />
  - [Isolating poly- and mono-cistronic reads](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-poly--and-mono-cistronic-reads) <br />
  - [Isolating reads with a SL](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-reads-with-a-sl) <br />
  - [Verify the SL is the same for all reads](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#verify-the-sl-is-the-same-for-all-reads) <br />
[Transcript evidence annotation and clean up](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#transcript-evidence-finding-and-clean-up) <br />

Software requirements: <br />
- emboss-6.6.0 fuzznuc <br />
- guppy-6.4.2 <br />
- manimap2.1 <br />
- meme-5.5.5 <br />
- ncbi-blast+-2.14.0 <br />
- samtools-1.20 <br />
- stringtie-2.2.1 <br />

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
input=path_to_fastq_to_analyze

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
Realign with [minimap2.sh](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/minimap.sh)<br />
and samtools (as previously described) <br />

### Verify the SL is the same for all reads
Isolating the 50 first nucleotides of each read <br />
```
seqkit subseq -r 1:50 monocistron_SL.fasta > 50first_nuc_SL.fasta
```
Using meme for motif search <br />
```
meme 50first_nuc_SL.fasta -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
```


## Transcript evidence finding and clean up

