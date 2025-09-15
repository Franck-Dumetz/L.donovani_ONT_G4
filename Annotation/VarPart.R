# --- Install required packages if needed ---
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

# --- Load libraries ---
library(car)
library(tidyverse)

# --- Load data ---
expr <- read.delim("PATH/transcript_counts_newgff.txt")
features <- read.delim("PATH/gene_features.tsv")

# --- Merge and prepare ---
dat <- merge(expr, features, by = "gene_id")

dat <- dat %>%
  mutate(
    PTU_ID = as.factor(PTU_ID),
    position_in_PTU = as.numeric(position_in_PTU),
    gene_length = as.numeric(gene_length),
    logCPM = log2(CPM + 1)  # log-transform CPM to stabilize variance
  ) %>%
  filter(!is.na(logCPM), is.finite(logCPM))

# --- Fit linear model with all predictors as fixed effects ---
fit_lm <- lm(logCPM ~ gene_length + position_in_PTU + PTU_ID, data = dat)

# --- Perform type III ANOVA ---
anova_table <- Anova(fit_lm, type = 3)

# --- Calculate partial eta-squared (variance explained) ---
ss_total <- sum(anova_table$`Sum Sq`)
eta_sq <- anova_table$`Sum Sq` / ss_total
names(eta_sq) <- rownames(anova_table)

# --- Print and visualize ---
print(round(eta_sq, 3))

# Clean up labels: replace underscores with spaces
clean_labels <- gsub("_", " ", names(eta_sq))

barplot(eta_sq,
        col = c("orangered", "steelblue", "red", "grey", "darkgreen"),
        ylab = "Proportion of Variance Explained",
        #main = "Fixed Effects Variance Partitioning",
        las = 1,                     # Rotate x labels vertically
        ylim = c(0, 1),              # Extend y-axis to full range (0 to 1)
        cex.names = 0.9,            # Slightly shrink font if needed
        names.arg = clean_labels)  # Ensure all labels are printed
