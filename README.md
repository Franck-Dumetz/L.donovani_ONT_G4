# _L. donovani_ UTR mapping (using ONT Direct RNA sequencing)

Table of content: <br />
[Basecalling ONT reads and mapping to the genome](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#basecalling-ont-reads-and-mapping-to-the-genome)<br />
[Isolating poly- and monocistrons with a splice leader](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#isolating-poly--and-mono-cistrons-with-a-splice-leader) <br />

Software requirements: <br />
• guppy-6.4.2 <br />
• manimap2.1 <br />
• ncbi-blast+-2.14.0 <br />

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

## Isolating poly- and mono-cistrons with a splice leader

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




