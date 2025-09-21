# Mapping G4 onto Ld1S genome and identifying rG4 in coding vs the UTRs of transcripts
Retreive G4 motifs from [canonical motifs described in "G-Quadruplex Identification in the Genome of Protozoan Parasites Points to Naphthalene Diimide Ligands as New Antiparasitic Agents"](https://bioinformatics.cruk.cam.ac.uk/G4Hunter/)![image](https://github.com/user-attachments/assets/721bd817-1a40-4355-a13d-cc6ce9e5283c) Bedrat, A., et al. (2016). "Re-evaluation of G-quadruplex propensity with G4Hunter." Nucleic Acids Res 44(4): 1746-1759. <br />
create a bed file for IGV visualisation and easy manipulation
```
for f in *.txt; do tail -n +2 "$f"; done > Ld1S_G4HunterSeeker.txt
```
```
awk 'NR>1 {print $1, $2-1, $3, "G4", $6, $5}' OFS="\t" Ld1S_G4HunterSeeker.txt > Ld1S_G4HunterSeeker.bed
```
Isolation of "CDS", "5'UTR" and "3'UTR" features from gff
```
input_gff=<path/to/gff>
grep -P '\tCDS\t' $input_gff > CDS.gff3
grep -P '\tfive_prime_UTR\t' $input_gff > five_prime_UTR.gff3
grep -P '\tthree_prime_UTR\t' $input_gff > three_prime_UTR.gff3
```
Convert to bed file
```
cut -f1,4,5,7 CDS.gff3 | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, ".", "0", $4}' > CDS.bed
cut -f1,4,5,7 five_prime_UTR.gff3 | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, ".", "0", $4}' > five_prime_UTR.bed
cut -f1,4,5,7 three_prime_UTR.gff3 | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, ".", "0", $4}' > three_prime_UTR.bed
```
Use BEDTools to find rG4 in the different sections of the transcript and get a count
```
input_bed=<path/to/bed>
bedtools intersect -a $input_bed -b CDS.bed -s -u | wc -l
bedtools intersect -a $input_bed -b five_prime_UTR.bed -s -u | wc -l
bedtools intersect -a $input_bed -b three_prime_UTR.bed -s -u | wc -l
```
## Gene expression
```
featureCounts -s 2 -T 8 -a /PATH/Ld1S_annotation_filtered_final_final.gtf -o Ld_3ONT_count.txt /PATH/Ld1S_3ONT_monocistron.bam
```
## Gene ontology analysis
Prepare the data
```
grep -A1 --no-group-separator -E '\.1\.p1([^0-9]|$)' linear_Ld1S_pep.fasta |sed 's/\*$//' > Ld1S_only.p1_linear.fasta
```
Use interproscan to find GO terms associated to each protein
```
module load interproscan
interproscan.sh -i AMA_pep_clean.fasta -f TSV -goterms -pa -cpu 16 -o AMA_blast_resultsGO.tsv
```
```
cut -f1,14 AMA_blast_resultsGO2.tsv | grep "GO:" AMA_spe_GO.txt
sed 's/|/,/g' AMA_spe_GO.txt | awk '{gsub(/\[[^]]*\]/, ""); gsub(/\{[^}]*\}/, ""); gsub(/\([^)]*\)/, ""); print}' > AMA_spe_GO_clean.txt
```
Use topGO in R like in [topGo](https://github.com/Franck-Dumetz/Ldonovani_UTR_mapping/blob/main/topGO.R) <br />

