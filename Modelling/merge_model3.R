## ==========================================
## Merge TE family VST with full merged data
## Files:
##  - BEAT_BM_c2_TEF_vst_9families.rds   (this is `tef`)
##  - BEAT_BM_c2_merged_full.csv         (clinical + scores + immune + secretome)
## Output:
##  - BEAT_BM_c2_merged_full_with_TEF.csv / .rds
## ==========================================

library(dplyr)
library(readr)
library(tibble)

## 1. Load TE family VST object -------------------------------------

tef <- readRDS("BEAT_BM_c2_TEF_vst_9families.rds")

# Make sure it's a data.frame
tef <- as.data.frame(tef)

cat("TEF dimensions:", dim(tef), "\n")
print(colnames(tef))

## Ensure we have a 'sample' column -------------------------------
## If it doesn't exist, create it from rownames
if (!"sample" %in% colnames(tef)) {
  message("No 'sample' column found in TEF object – using rownames as sample IDs.")
  tef <- tef %>%
    rownames_to_column(var = "sample")
}

# Convert TE family columns to numeric (leave 'sample' as character)
tef <- tef %>%
  mutate(
    sample = as.character(sample),
    across(-sample, as.numeric)
  )

# Quick peek
tef[1:5, ]

## 2. Load main merged dataset --------------------------------------

full <- read_csv("BEAT_BM_c2_merged_full.csv", show_col_types = FALSE)

cat("Full merged dimensions:", dim(full), "\n")
head(full[, 1:5])

# Key column should be: dbgap_rnaseq_sample
if (!"dbgap_rnaseq_sample" %in% colnames(full)) {
  stop("Column 'dbgap_rnaseq_sample' not found in BEAT_BM_c2_merged_full.csv")
}

## 3. Prepare TEF table for join ------------------------------------

# Rename 'sample' -> 'dbgap_rnaseq_sample' to match `full`
tef <- as.data.frame(tef)

tef$dbgap_rnaseq_sample <- tef$sample
tef$sample <- NULL


# Optional: check overlap
common_ids <- intersect(full$dbgap_rnaseq_sample, tef$dbgap_rnaseq_sample)
cat("Common samples between FULL and TEF:", length(common_ids), "\n")

## 4. Merge ----------------------------------------------------------

merged_with_tef <- full %>%
  left_join(tef, by = "dbgap_rnaseq_sample")

cat("Merged + TEF dimensions:", dim(merged_with_tef), "\n")
merged_with_tef[1:3, 1:10]

## 5. Save output ----------------------------------------------------

write_csv(merged_with_tef, "BEAT_BM_c2_merged_full_with_TEF.csv")
saveRDS(merged_with_tef, "BEAT_BM_c2_merged_full_with_TEF.rds")

cat("Saved:\n",
    " - BEAT_BM_c2_merged_full_with_TEF.csv\n",
    " - BEAT_BM_c2_merged_full_with_TEF.rds\n")
