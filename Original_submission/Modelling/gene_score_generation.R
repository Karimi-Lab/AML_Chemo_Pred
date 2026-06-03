## ============================================================
## BEAT_BM_c2: Generate singscore from TPM
## Inputs:
##   - BEAT_BM_c2_ge_tpm.txt        (genes x samples, TPM)
##   - Cell-type-composition.csv    (columns = gene sets)
##   - Immune_profile_gene_set.csv  (columns = gene sets)
##   - Ageing.csv                   (columns = gene sets)
## Output:
##   - BEAT_BM_c2_singscore_fromTPM.rds (samples x signatures)
##   - BEAT_BM_c2_singscore_fromTPM.csv
## ============================================================

library(dplyr)
library(readr)
library(singscore)

## 0. Config --------------------------------------------------
expr_file   <- "BEAT_BM_c2_ge_tpm.txt"
out_prefix  <- "BEAT_BM_c2"

celltype_file <- "Cell-type-composition.csv"
immune_file   <- "Immune_profile_gene_set.csv"
ageing_file   <- "Ageing.csv"

## 1. Load TPM expression matrix (genes x samples) ------------
expr_tpm <- read.delim(
  expr_file,
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

cat("TPM expression matrix dimensions:", dim(expr_tpm), "\n")
print(expr_tpm[1:5, 1:5])

# Convert to numeric matrix and standardise gene names to uppercase
expr_tpm <- as.matrix(expr_tpm)
rownames(expr_tpm) <- toupper(rownames(expr_tpm))

## 2. Helper: read gene sets from a CSV (columns = gene sets) ---------
## Each column = one geneset; entries are gene symbols (with NAs).

read_genesets_from_csv <- function(file) {
  df <- read.csv(file,
                 header           = TRUE,
                 stringsAsFactors = FALSE,
                 check.names      = FALSE)
  
  gs_list <- lapply(df, function(col) {
    genes <- unique(na.omit(col))
    genes <- genes[genes != ""]
    toupper(genes)
  })
  
  names(gs_list) <- colnames(df)
  gs_list
}

## 3. Load gene sets from the three files -----------------------------
gs_celltype <- read_genesets_from_csv(celltype_file)
gs_immune   <- read_genesets_from_csv(immune_file)
gs_ageing   <- read_genesets_from_csv(ageing_file)

cat("Number of cell-type gene sets:", length(gs_celltype), "\n")
cat("Number of immune gene sets:",    length(gs_immune),   "\n")
cat("Number of ageing gene sets:",    length(gs_ageing),   "\n")

# Combine into one big list of gene sets
gene_sets <- c(gs_celltype, gs_immune, gs_ageing)
cat("Total number of gene sets:", length(gene_sets), "\n")
cat("Gene set names:\n")
print(names(gene_sets))

## 4. Rank genes (required input for singscore) -----------------------
# expr_tpm: genes x samples
rank_data <- singscore::rankGenes(expr_tpm)
cat("rank_data class:", class(rank_data), "\n")
cat("rank_data dimensions:", dim(rank_data), "\n")

## 5. Compute singscore for each gene set -----------------------------
# We treat all gene sets as "up" signatures (no downSet).

scores_mat <- sapply(names(gene_sets), function(gs_name) {
  gs_genes <- gene_sets[[gs_name]]
  
  # overlap with genes present in the expression matrix
  common <- intersect(gs_genes, rownames(rank_data))
  
  if (length(common) < 5) {
    message("[", out_prefix, "] Gene set ", gs_name, ": only ",
            length(common), " genes in common with expression; setting scores to NA.")
    return(rep(NA_real_, ncol(rank_data)))
  }
  
  sc <- singscore::simpleScore(
    rankData    = rank_data,  # MUST be RankData, not a matrix
    upSet       = common
  )
  
  # TotalScore is a numeric vector of length = number of samples
  sc$TotalScore
})

cat("scores_mat dimensions (samples x signatures?):", dim(scores_mat), "\n")

## 6. Convert scores matrix to a clean data.frame ---------------------
# sapply above returns a matrix with rows = samples, cols = gene sets
scores_df <- as.data.frame(scores_mat, stringsAsFactors = FALSE)

# Add sample column from the column names of rank_data (same as expr_tpm)
scores_df$sample <- colnames(rank_data)

# Move 'sample' to the first column
scores_df <- scores_df %>%
  relocate(sample)

# Clean column names (valid for modelling)
colnames(scores_df)[-1] <- make.names(colnames(scores_df)[-1], unique = TRUE)

cat("Final singscore data.frame dimensions (samples x features):",
    dim(scores_df), "\n")
print(scores_df[1:5, 1:5])

## 7. Save results -----------------------------------------------------
rds_out  <- paste0(out_prefix, "_singscore_fromTPM.rds")
csv_out  <- paste0(out_prefix, "_singscore_fromTPM.csv")

saveRDS(scores_df, file = rds_out)
write.csv(scores_df, file = csv_out, row.names = FALSE)

cat("Singscore generation complete. Files written:\n",
    "  - ", rds_out, "\n",
    "  - ", csv_out, "\n", sep = "")
