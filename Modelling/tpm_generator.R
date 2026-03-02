## ============================================================
## BEAT_BM_c2: Generate TPM for all genes from raw counts
## Inputs:
##   - BEAT_BM_c2_gene_counts.txt   (genes x samples, raw counts)
##   - gene_length_kb.Rda           (named numeric, exonic length in kb)
## Outputs:
##   - BEAT_BM_c2_ge_tpm.txt        (genes x samples, TPM)
##   - BEAT_BM_c2_ge_log2tpm.txt    (genes x samples, log2(TPM+1))
##   - BEAT_BM_c2_ge_tpm.rds
##   - BEAT_BM_c2_ge_log2tpm.rds
## ============================================================

library(dplyr)
library(tibble)

## 1. Load raw gene counts (genes x samples) -------------------------

counts_c2 <- read.delim(
  "/Users/mehdi/Library/CloudStorage/OneDrive-King'sCollegeLondon(2)/Staffs/Sila_Gerlevic/AML/data/New_data/BEAT_BM_c2_gene_counts.txt",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

cat("Raw counts matrix dimensions:", dim(counts_c2), "\n")
print(counts_c2[1:5, 1:5])


## 2. Normalise gene names (uppercase) and collapse duplicates -------

# Move rownames to a column, uppercase, and sum duplicates
counts_c2_clean <- counts_c2 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  as.data.frame()

# Restore as matrix with unique uppercase rownames
rownames(counts_c2_clean) <- counts_c2_clean$gene
counts_c2_clean$gene <- NULL

cat("Counts matrix after collapsing duplicate genes:", dim(counts_c2_clean), "\n")
print(counts_c2_clean[1:5, 1:5])

# Sanity check: no duplicate rownames
stopifnot(!any(duplicated(rownames(counts_c2_clean))))


## 3. Load gene lengths in kb (named numeric vector) -----------------

gene_length_kb <- readRDS("gene_length_kb.Rda")
cat("Example of gene_length_kb:\n")
print(head(gene_length_kb))

# Uppercase gene names to match counts
names(gene_length_kb) <- toupper(names(gene_length_kb))

# Intersect genes present in both counts and length vector
common_genes <- intersect(rownames(counts_c2_clean), names(gene_length_kb))
cat("Number of genes with both counts and length:", length(common_genes), "\n")

# Subset counts and lengths
counts_c2_common      <- counts_c2_clean[common_genes, , drop = FALSE]
gene_length_kb_common <- gene_length_kb[common_genes]


## 4. Helper function: counts -> TPM using length in kb --------------

counts_to_tpm <- function(count_mat, gene_length_kb) {
  # count_mat: genes x samples (numeric matrix/data.frame)
  # gene_length_kb: named numeric vector, length in kb, names = genes
  
  # Ensure numeric matrix
  count_mat <- as.matrix(count_mat)
  
  # Reorder length vector to match gene order
  gene_length_kb <- gene_length_kb[rownames(count_mat)]
  
  # rate = counts / length_kb
  rate <- sweep(count_mat, 1, gene_length_kb, FUN = "/")
  
  # Sum over genes per sample
  scale <- colSums(rate, na.rm = TRUE)
  
  # TPM = rate / sum(rate) * 1e6
  tpm <- sweep(rate, 2, scale, FUN = "/") * 1e6
  
  return(tpm)
}


## 5. Compute TPM for ALL common genes -------------------------------

tpm_c2 <- counts_to_tpm(counts_c2_common, gene_length_kb_common)

cat("TPM matrix dimensions:", dim(tpm_c2), "\n")
print(tpm_c2[1:5, 1:5])


## 6. Optional: log2(TPM + 1) ----------------------------------------

log2_tpm_c2 <- log2(tpm_c2 + 1)

cat("log2(TPM+1) example:\n")
print(log2_tpm_c2[1:5, 1:5])


## 7. Save TPM matrices ----------------------------------------------

# Save as tab-delimited text (genes x samples)
write.table(
  tpm_c2,
  file      = "BEAT_BM_c2_ge_tpm.txt",
  quote     = FALSE,
  sep       = "\t",
  row.names = TRUE,
  col.names = NA
)

write.table(
  log2_tpm_c2,
  file      = "BEAT_BM_c2_ge_log2tpm.txt",
  quote     = FALSE,
  sep       = "\t",
  row.names = TRUE,
  col.names = NA
)

# Save as R objects
saveRDS(tpm_c2,      file = "BEAT_BM_c2_ge_tpm.rds")
saveRDS(log2_tpm_c2, file = "BEAT_BM_c2_ge_log2tpm.rds")

cat("TPM generation complete. Files written:\n",
    "  - BEAT_BM_c2_ge_tpm.txt / .rds\n",
    "  - BEAT_BM_c2_ge_log2tpm.txt / .rds\n")



library(dplyr)

## 1. Secretome gene list -----------------------------------------

secretome_genes <- c(
  "CCL3","CX3CL1","HGF","FGF4","CXCL2","MMP2","IL11","MMP9","CST3","TGFB1",
  "SERPINE1","CXCL12","CCL2","IL2","IFNG","VEGFA","HBEGF","IL4","IL5","FGF1",
  "IL1A","SERPINc2","ARG1","TEK","IL1B","CCL4","IL6","IL10","IL1RN","IL33",
  "FGF2","CXCL9","EGF","MMP3","IL18","CSTB","CCL5","TGFA","CXCR1","CXCL5",
  "CXCL1","IL15","CSF2","IL25","IL12A","IL13","CXCL10","CXCL11","IL8","IFNB1",
  "BSG","CXCR2","FGF3","MMP1","ELANE","CFD","IFNA1","TNF"
)

## 2. Read TPM matrix (genes x samples) ----------------------------

tpm_c2 <- read.delim(
  "BEAT_BM_c2_ge_tpm.txt",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

cat("TPM matrix dimensions:", dim(tpm_c2), "\n")
tpm_c2[1:5, 1:5]

## 3. Harmonise gene names and select secretome genes --------------

# TPM file from previous step should already have uppercase rownames,
# but we enforce it just to be safe
rownames(tpm_c2) <- toupper(rownames(tpm_c2))
secretome_genes_up <- toupper(secretome_genes)

sec_genes_in      <- intersect(secretome_genes_up, rownames(tpm_c2))
sec_genes_missing <- setdiff(secretome_genes_up, rownames(tpm_c2))

cat("Secretome genes found in TPM matrix:",
    length(sec_genes_in), "/", length(secretome_genes_up), "\n")

if (length(sec_genes_missing) > 0) {
  cat("Missing secretome genes (not in TPM matrix):\n",
      paste(sec_genes_missing, collapse = ", "), "\n")
}

# Subset TPM to secretome genes
sec_tpm_c2 <- tpm_c2[sec_genes_in, , drop = FALSE]   # genes x samples

## 4. Convert to samples x genes and tidy ---------------------------

secretome_tpm_c2 <- as.data.frame(t(sec_tpm_c2))
secretome_tpm_c2$sample <- rownames(secretome_tpm_c2)
rownames(secretome_tpm_c2) <- NULL

# Make gene column names syntactically valid for modelling
names(secretome_tpm_c2) <- make.names(names(secretome_tpm_c2), unique = TRUE)

cat("Secretome TPM matrix (samples x genes):", dim(secretome_tpm_c2), "\n")
secretome_tpm_c2[1:5, 1:5]

## 5. Save secretome TPM --------------------------------------------

saveRDS(secretome_tpm_c2, "BEAT_BM_c2_secretome_tpm.rds")
write.csv(secretome_tpm_c2, "BEAT_BM_c2_secretome_tpm.csv", row.names = FALSE)

cat("Secretome TPM saved as:\n",
    "  - BEAT_BM_c2_secretome_tpm.rds\n",
    "  - BEAT_BM_c2_secretome_tpm.csv\n")

