## ============================================================
## BEAT_BM_c2: VST for gene counts + TE family counts
## ============================================================

library(DESeq2)
library(dplyr)
library(tibble)

## -----------------------------
## 1. Load gene count matrix
## -----------------------------
gene_counts <- read.delim(
  "/Users/mehdi/Library/CloudStorage/OneDrive-King'sCollegeLondon(2)/Staffs/Sila_Gerlevic/AML/data/New_data/BEAT_BM_c2_gene_counts.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

message("Gene counts: ", nrow(gene_counts), " genes x ", ncol(gene_counts), " samples")

## -----------------------------
## 2. Load TE family counts Rda
## -----------------------------
te_env <- new.env()
load("/Users/mehdi/Library/CloudStorage/OneDrive-King'sCollegeLondon(2)/Staffs/Sila_Gerlevic/AML/data/New_data/BEAT_BM_TE_family_counts.Rda", envir = te_env)

# Pick first matrix-like object
te_obj <- NULL
for (nm in ls(te_env)) {
  obj <- get(nm, envir = te_env)
  if (is.matrix(obj) || is.data.frame(obj)) {
    te_obj <- obj
    message("Using TE object: ", nm)
    break
  }
}
stopifnot(!is.null(te_obj))

te_df <- as.data.frame(te_obj, check.names = FALSE)

## ---- FIX: assign TE_family as rownames ----
if (!"TE_family" %in% colnames(te_df)) {
  stop("TE_family column missing in Rda file.")
}

te_fam_names <- te_df$TE_family
te_mat <- te_df[, colnames(te_df) != "TE_family", drop = FALSE]

# Convert to numeric
te_mat <- as.data.frame(lapply(te_mat, function(x) as.numeric(trimws(x))),
                        check.names = FALSE)
te_counts <- as.matrix(te_mat)

# Assign corrected rownames
rownames(te_counts) <- te_fam_names
rownames(te_counts) <- paste0("TEF_", rownames(te_counts))

message("TE families: ", nrow(te_counts), " families")

print(te_counts[1:5, 1:5])

## -----------------------------
## 3. Align samples
## -----------------------------
common_samples <- intersect(colnames(gene_counts), colnames(te_counts))

message("Common samples: ", length(common_samples))

gene_counts <- gene_counts[, common_samples, drop = FALSE]
te_counts   <- te_counts[, common_samples, drop = FALSE]

## -----------------------------
## 4. Combine matrices
## -----------------------------
all_counts <- rbind(gene_counts, te_counts)

## -----------------------------
## 5. Run VST
## -----------------------------
# DESeq2 requires a metadata frame, even if unused
coldata <- data.frame(
  dummy = rep(1, ncol(all_counts))  # one row per sample
)
rownames(coldata) <- colnames(all_counts)

dds <- DESeqDataSetFromMatrix(
  countData = round(all_counts),
  colData   = coldata,
  design    = ~ 1   # no group comparison
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)


## 8. Split back into gene vs TE family VST (optional) -------------
is_te <- grepl("^TEF_", rownames(vsd_mat))

vsd_genes <- vsd_mat[!is_te, , drop = FALSE]
vsd_te    <- vsd_mat[ is_te, , drop = FALSE]

cat("VST genes:", dim(vsd_genes), "\n")
cat("VST TE families:", dim(vsd_te), "\n")

## 9. Save objects --------------------------------------------------
saveRDS(vsd_mat,  "BEAT_BM_c2_gene_TEF_vst_all.rds")
saveRDS(vsd_genes,"BEAT_BM_c2_gene_vst.rds")
saveRDS(vsd_te,   "BEAT_BM_c2_TEF_vst.rds")





## ============================================
## Make filtered TE-family VST file (9 families)
## ============================================

library(dplyr)

## 0) Load TE-family VST matrix -------------------------------------
## vsd_te: rows = TE families (e.g. "TEF_LINE:L1"), cols = samples

# If vsd_te is already in your session, skip this.
# Otherwise, uncomment and adjust the path:
# vsd_te <- readRDS("BEAT_BM_c2_TEF_vst_all.rds")

stopifnot(exists("vsd_te"))

cat("Full TE VST matrix dimensions (families x samples):\n")
print(dim(vsd_te))
cat("Example rows:\n")
print(head(rownames(vsd_te), 10))

## 1) TE families to keep (exact names as in assay) -----------------

te_keep_prefixed <- c(
  "TEF_LINE:L2",
  "TEF_LTR:ERVL",
  "TEF_SINE:MIR",
  "TEF_LTR:ERV1",
  "TEF_LTR:ERVL-MaLR",
  "TEF_LINE:L1",
  "TEF_LTR:ERVK",
  "TEF_SINE:Alu",
  "TEF_LINE:CR1"
)

## 2) Subset vsd_te to these families -------------------------------

te_present <- intersect(te_keep_prefixed, rownames(vsd_te))
te_missing <- setdiff(te_keep_prefixed, rownames(vsd_te))

cat("\nRequested TE families:", length(te_keep_prefixed), "\n")
cat("Found in vsd_te:        ", length(te_present), "\n")

if (length(te_missing) > 0) {
  cat("Missing TE families (not in vsd_te):\n",
      paste(te_missing, collapse = ", "), "\n")
}

vsd_te_sub <- vsd_te[te_present, , drop = FALSE]

cat("\nSubmatrix dimensions (selected TE families x samples):\n")
print(dim(vsd_te_sub))

## 3) Remove 'TEF_' prefix in the output ---------------------------

rownames(vsd_te_sub) <- sub("^TEF_", "", rownames(vsd_te_sub))

## 4) Transpose to samples x TE-families for merging ----------------

te_vst_filtered <- t(vsd_te_sub) %>% as.data.frame()

# Add sample column for joins later
te_vst_filtered$sample <- rownames(te_vst_filtered)
rownames(te_vst_filtered) <- NULL

# Make column names syntactically valid (just in case),
# but the visible labels will still be like "LINE:L1", "LTR:ERVK", etc.
names(te_vst_filtered) <- make.names(names(te_vst_filtered), unique = TRUE)

cat("\nFinal filtered TE VST matrix (samples x TE variables):\n")
print(dim(te_vst_filtered))
cat("Column names:\n")
print(names(te_vst_filtered))

## 5) Write filtered files ------------------------------------------

saveRDS(te_vst_filtered, "BEAT_BM_c2_TEF_vst_9families.rds")
write.csv(te_vst_filtered,
          "BEAT_BM_c2_TEF_vst_9families.csv",
          row.names = FALSE)

cat("\nWritten files:\n",
    "  - BEAT_BM_c2_TEF_vst_9families.rds\n",
    "  - BEAT_BM_c2_TEF_vst_9families.csv\n")




