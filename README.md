# Ld1S_genome
Making genome assembly of Ld1S2D and further comparison with BPK282A2 (TritrypDB v63)

#SNP density

* MUMmer SNP comparition
Using nucmer to identify SNP and position

/usr/local/packages/mummer-3.23/nucmer -p Ld1Svs282_WG_nucmer --mum <PATH-TO-REF-GENOME>/TriTrypDB-63_LdonovaniBPK282A1_Genome.fasta <PATH-TO-QUERY_GENOME>/Ld1S_assembly_final.fasta

/usr/local/packages/mummer-3.23/delta-filter -r -q Ld1Svs282_WG_nucmer.delta > Ld1Svs282_WG_nucmer.filter

/usr/local/packages/mummer-3.23/show-snps -Clr Ld1Svs282_WG_nucmer.filter > Ld1Svs282_WG_nucmer.snps

Filtering Ld1Svs282_WG_nucmer.snps to extract only SNPs (no insertion, no deletion)

tail -n +6 Ld1Svs282_WG_nucmer.snps | awk '{if($2!=".") print}' | awk '{if($3!=".") print $1"\t"$14}' > SNP_density.tsv

Importing the dataset into R to plot the density for all chromosomes

#!/usr/bin/env Rscript

library(ggplot2)

colnames(SNP_density) <- c("Position", "Chr")
SNP_density$Position <- as.numeric(SNP_density$Position)

head(SNP_density)

ggplot(aes(Position, group = Chr), data=SNP_density) +
    geom_histogram(binwidth = 5000) + 
    #facet_wrap(~Chr, scales = "free") +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
