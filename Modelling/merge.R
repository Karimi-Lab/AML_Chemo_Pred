## ================================================================
## MERGE PIPELINE – BEAT_BM C1
## Clinical + Singscore + CIBERSORT LM22 + Secretome TPM
## Output: BEAT_BM_c2_merged_full.csv / .rds
## ================================================================

library(dplyr)
library(readr)

## ---------------------------------------------------------------
## 1. Load clinical data
## ---------------------------------------------------------------

clin_c2 <- read_csv("aml_baseline_clinical_beats_BM_c2_cyto.csv")

cat("Clinical rows:", nrow(clin_c2), "\n")
clin_c2[1:3, 1:5]

## This must match the sample identifiers used in expression data:
## The merge key is: dbgap_rnaseq_sample


## ---------------------------------------------------------------
## 2. Load Singscore output (samples x signatures)
## ---------------------------------------------------------------

singscore_c2 <- read_csv("BEAT_BM_c2_singscore_fromTPM.csv")

# Ensure sample column is named "sample"
names(singscore_c2)[1] <- "sample"

cat("Singscore rows:", nrow(singscore_c2), "\n")
singscore_c2[1:3, 1:5]


## ---------------------------------------------------------------
## 3. Load CIBERSORT LM22 results
## ---------------------------------------------------------------

ciber_c2 <- read_csv("BEAT_BM_c2_CIBERSORT_LM22_results.csv")

# Usually CIBERSORT writes 'Mixture' or the sample ID in first column
names(ciber_c2)[1] <- "sample"

cat("CIBERSORT rows:", nrow(ciber_c2), "\n")
ciber_c2[1:3, 1:5]


## ---------------------------------------------------------------
## 4. Load Secretome TPM values
## ---------------------------------------------------------------

secretome_c2 <- read_csv("BEAT_BM_c2_secretome_tpm.csv")

# Ensure sample ID column is named "sample"
names(secretome_c2)[ncol(secretome_c2)] <- "sample"

cat("Secretome rows:", nrow(secretome_c2), "\n")
secretome_c2[1:3, 1:5]


## ---------------------------------------------------------------
## 5. Perform the merges
## ---------------------------------------------------------------

merged_c2 <- clin_c2 %>%
  left_join(singscore_c2, by = c("dbgap_rnaseq_sample" = "sample")) %>%
  left_join(ciber_c2,     by = c("dbgap_rnaseq_sample" = "sample")) %>%
  left_join(secretome_c2, by = c("dbgap_rnaseq_sample" = "sample"))

cat("\nMerged dimensions:", dim(merged_c2), "\n")
merged_c2[1:3, 1:10]


## ---------------------------------------------------------------
## 6. Save output
## ---------------------------------------------------------------

write_csv(merged_c2, "BEAT_BM_c2_merged_full.csv")
saveRDS(merged_c2, "BEAT_BM_c2_merged_full.rds")

cat("\nSaved merged dataset as:\n",
    " - BEAT_BM_c2_merged_full.csv\n",
    " - BEAT_BM_c2_merged_full.rds\n")
