# Install if needed
# install.packages(c("readr", "dplyr", "pheatmap", "Hmisc"))

library(readr)
library(dplyr)
library(pheatmap)
library(Hmisc)

# 1. Load data
df <- read_csv("BEAT_BM_c1_merged_full_with_TEF.csv")

# 2. Select ERV features
erv_features <- c(
  "LTR.ERV1",
  "LTR.ERVL",
  "LTR.ERVL.MaLR",
  "LTR.ERVK"
)

# 3. Select transcriptomic/program features
transcriptomic_features <- c(
  "HSC.like",
  "HSC.Prog",
  "Progenitor.like",
  "IFNA1",
  "IFNB1",
  "CXCL10",
  "CXCL11",
  "TNF",
  "IL10",
  "SASP",
  "Cellular_senescence",
  "T.cell.senescence",
  "Proliferation",
  "Exhausted.T.cells",
  "Cytotoxicity"
)

# 4. Keep available columns only
erv_features <- intersect(erv_features, colnames(df))
transcriptomic_features <- intersect(transcriptomic_features, colnames(df))

# 5. Convert to numeric
cor_df <- df %>%
  select(all_of(c(erv_features, transcriptomic_features))) %>%
  mutate(across(everything(), as.numeric))

# 6. Spearman correlation: ERVs vs transcriptomic features
cor_results <- rcorr(
  as.matrix(cor_df[, c(erv_features, transcriptomic_features)]),
  type = "spearman"
)

cor_mat <- cor_results$r[erv_features, transcriptomic_features]
p_mat <- cor_results$P[erv_features, transcriptomic_features]

# 7. Save outputs
write.csv(cor_mat, "ERV_transcriptomic_spearman_correlation_matrix.csv")
write.csv(p_mat, "ERV_transcriptomic_spearman_pvalue_matrix.csv")

# 8. Plot heatmap
pdf("ERV_transcriptomic_correlation_heatmap.pdf", width = 10, height = 5)
pheatmap(
  cor_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Spearman correlation: ERV features vs transcriptomic programmes"
)
dev.off()