# Ld1S_genome
Making genome assembly of Ld1S2D and further comparison with BPK282A2 (TritrypDB v63)


* Genome assembly strategies
   ** Using Flye
export PATH=/usr/local/packages/flye-2.9/bin:$PATH  
/usr/local/packages/flye-2.9/bin/flye --pacbio-hifi Ld1S2D_DNA_filtered.ccs.fastq.gz --genome-size 33m --out-dir /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S2D_genome/Flye -t 16

   ** Using Canu


* Checking genome assembly quality compared to LdBPK282A2
   ** using MUMmer

/usr/local/packages/mummer-3.23/nucmer -p Ld1Svs282_WG_nucmer --mum <PATH-TO-REF-GENOME>/TriTrypDB-63_LdonovaniBPK282A1_Genome.fasta <PATH-TO-QUERY_GENOME>/Ld1S_assembly_final.fasta
/usr/local/packages/mummer-3.23/mummerplot --png --filter --color --layout --prefix=Ld1Svs282_WG_nucmer Ld1Svs282_WG_nucmer.delta -R <PATH-TO-REF-GENOME>/TriTrypDB-63_LdonovaniBPK282A1_Genome.fasta -Q <PATH-TO-QUERY_GENOME>/Ld1S_assembly_final.fasta

* Rearrangement
/usr/local/packages/mummer-3.23/show-diff Ld1Svs282_WG_nucmer.delta > Ld1Svs282_WG_nucmer.rearrangement

* Repeats

/usr/local/packages/mummer-3.23/nucmer --maxmatch --nosimplify --prefix=Ld1Svs282_WG_nucmermax /local/projects-t3/SerreDLab-3/fdumetz/Genomes/LdBPK282A2/TriTrypDB-63_LdonovaniBPK282A1_Genome.fasta /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_genome/Flye_scaffold/Ld1S_assembly_final.fasta
/usr/local/packages/mummer-3.23/show-coords -r Ld1Svs282_WG_nucmermax.delta > Ld1Svs282_WG_nucmermax.coords

* SNP density

    ** MUMmer SNP comparition
Using nucmer to identify SNP and position

/usr/local/packages/mummer-3.23/nucmer -p Ld1Svs282_WG_nucmer --mum <PATH-TO-REF-GENOME>/TriTrypDB-63_LdonovaniBPK282A1_Genome.fasta <PATH-TO-QUERY_GENOME>/Ld1S_assembly_final.fasta

/usr/local/packages/mummer-3.23/delta-filter -r -q Ld1Svs282_WG_nucmer.delta > Ld1Svs282_WG_nucmer.filter

/usr/local/packages/mummer-3.23/show-snps -Clr Ld1Svs282_WG_nucmer.filter > Ld1Svs282_WG_nucmer.snps

   ** Filtering Ld1Svs282_WG_nucmer.snps to extract only SNPs (no insertion, no deletion)

tail -n +6 Ld1Svs282_WG_nucmer.snps | awk '{if($2!=".") print}' | awk '{if($3!=".") print $1"\t"$14}' > SNP_density.tsv

   ** Importing the dataset into R to plot the density for all chromosomes

#!/usr/bin/env Rscript

library(ggplot2)

colnames(SNP_density) <- c("Position", "Chr")
SNP_density$Position <- as.numeric(SNP_density$Position)

head(SNP_density)

ggplot(aes(Position, group = Chr), data=SNP_density) +
    geom_histogram(binwidth = 5000) + 
    #facet_wrap(~Chr, scales = "free") +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE))


----------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------
* tRNA destection using tRNA scan

/usr/local/packages/trnascan-se-2.0.3/bin/eufindtRNA -r <PATH>/Ld1S_assembly_final.fasta > Ld1S_tRNA_strict.csv 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------

* snoRNA detection using snoReport2.0

export SNOREPORTMODELS="/usr/local/packages/snoreport-2.0/models"   #### make  sure the "models" folder is readable or it won't work
/usr/local/packages/snoreport-2.0/snoreport_2 -i /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_genome/Flye_scaffold/Ld1S_assembly_final.fasta -CD -HACA -o Ld1S_snoRNA --PS

----------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------
* Using BUSCO

export PATH="/usr/local/packages/augustus-3.4.0/bin:$PATH"
export PATH="/usr/local/packages/augustus-3.4.0/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/usr/local/packages/augustus-3.4.0/configs/"
export PATH="/usr/local/packages/metaeuk-6-a5d39d9/bin:$PATH"
busco -m genome -i <PATH-TO-QUERY_GENOME>/Ld1S_assembly_final.fasta --auto-lineage-euk --long -o busco  -f





