# Load required libraries
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("topGO")
library(topGO)

# ---- Step 1: Load Gene Lists ----
# Read gene lists for two datasets
dataset1_genes <- read.table("/Volumes/projects-t3/SerreDLab-3/fdumetz/Leishmania/GO_analysis/AMA_spe_GO_clean2.txt", stringsAsFactors = FALSE)$V1
dataset2_genes <- read.table("/Volumes/projects-t3/SerreDLab-3/fdumetz/Leishmania/GO_analysis/PRO_spe_GO_clean2.txt", stringsAsFactors = FALSE)$V1

# Read background gene list (all detected genes in the experiment)
all_genes <- read.table("/Volumes/projects-t3/SerreDLab-3/fdumetz/Leishmania/GO_analysis/All_GO_clean2.txt", stringsAsFactors = FALSE)$V1

# ---- Step 2: Create Gene Factor Lists for `topGO` ----
# Initialize all genes as not significant (0), mark dataset genes as significant (1)
geneList1 <- factor(as.integer(all_genes %in% dataset1_genes))
names(geneList1) <- all_genes

geneList2 <- factor(as.integer(all_genes %in% dataset2_genes))
names(geneList2) <- all_genes

# Load gene-to-GO mappings (Format: gene_ID -> GO terms)
gene2GO <- readMappings(file = "/Volumes/projects-t3/SerreDLab-3/fdumetz/Leishmania/GO_analysis/onlyp1_onlyGO_clean2.txt")

# Create `topGOdata` objects for biological process (BP) ontology
GOdata1 <- new("topGOdata", ontology = "BP", allGenes = geneList1, annot = annFUN.gene2GO, gene2GO = gene2GO)
GOdata2 <- new("topGOdata", ontology = "BP", allGenes = geneList2, annot = annFUN.gene2GO, gene2GO = gene2GO)

# ---- Step 43: Perform GO Enrichment Analysis ----
# Run Fisher's exact test for enrichment
result1 <- runTest(GOdata1, algorithm = "weight01", statistic = "Fisher")
result2 <- runTest(GOdata2, algorithm = "weight01", statistic = "Fisher")

# Extract significant GO terms (adjust topNodes as needed)
table1 <- GenTable(GOdata1, result1, topNodes = 20)
table2 <- GenTable(GOdata2, result2, topNodes = 20)

# Save enrichment results
write.table(table1, "GO_enrichment_dataset1.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(table2, "GO_enrichment_dataset2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Step 4.1: Identify Shared & Unique GO Terms ----
shared_GO <- intersect(table1$GO.ID, table2$GO.ID)
unique_GO1 <- setdiff(table1$GO.ID, table2$GO.ID)
unique_GO2 <- setdiff(table2$GO.ID, table1$GO.ID)

# Save results
write.table(shared_GO, "shared_GO_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = "GO.ID")
write.table(unique_GO1, "unique_GO_dataset1.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = "GO.ID")
write.table(unique_GO2, "unique_GO_dataset2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = "GO.ID")

# Print summary
cat("GO Enrichment Analysis Completed\n")
cat(length(shared_GO), "shared GO terms found.\n")
cat(length(unique_GO1), "unique GO terms in Dataset 1.\n")
cat(length(unique_GO2), "unique GO terms in Dataset 2.\n")
