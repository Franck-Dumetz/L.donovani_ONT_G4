# _L. donovani_ UTR mapping (using ONT Direct RNA sequencing)

Table of content: <br />
[Basecalling ONT reads and mapping to the genome](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/README.md#basecalling-ont-reads-and-mapping-to-the-genome)<br />
[Isolating poly- and monocistrons with a splice leader]() <br />

Software requirements: <br />
• guppy-6.4.2 <br />
• manimap2.1 <br />

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

## Isolating poly- and monocistrons with a splice leader







