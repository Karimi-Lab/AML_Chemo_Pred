## ============================================================
## Libraries
## ============================================================
library(dplyr)
library(stringr)

## ============================================================
## Function to add cytogenetic flags and drop karyotype
## ============================================================

add_cytogenetic_flags <- function(df) {
  
  if (!"karyotype" %in% colnames(df)) {
    stop("Column 'karyotype' not found in the data frame.")
  }
  
  df %>%
    # keep raw text if you want it later – comment this line out if not needed
    rename(karyotype_raw = karyotype) %>%
    mutate(
      karyotype_clean = if_else(is.na(karyotype_raw), "", karyotype_raw),
      
      ## Core favourable groups
      cyto_t_8_21  = str_detect(karyotype_clean, "t\\(8;21\\)"),
      cyto_inv16   = str_detect(karyotype_clean, "inv\\(16\\)") |
        str_detect(karyotype_clean, "t\\(16;16\\)"),
      cyto_cbf     = cyto_t_8_21 | cyto_inv16,
      
      ## APL
      cyto_t_15_17 = str_detect(karyotype_clean, "t\\(15;17\\)"),
      
      ## Selected adverse abnormalities
      cyto_inv3    = str_detect(karyotype_clean, "inv\\(3\\)") |
        str_detect(karyotype_clean, "t\\(3;3\\)"),
      cyto_t_6_9   = str_detect(karyotype_clean, "t\\(6;9\\)"),
      cyto_5_7_abn = str_detect(karyotype_clean, "-5")  |
        str_detect(karyotype_clean, "del\\(5") |
        str_detect(karyotype_clean, "-7")  |
        str_detect(karyotype_clean, "del\\(7"),
      cyto_plus8   = str_detect(karyotype_clean, "\\+8"),
      
      ## Normal vs abnormal
      cyto_normal  = karyotype_clean %in% c("46,XY[20]", "46,XX[20]"),
      
      ## Very rough "complex" definition: lots of commas/semicolons, not normal
      cyto_complex = !cyto_normal &
        (str_count(karyotype_clean, ",") + str_count(karyotype_clean, ";") >= 5)
    ) %>%
    # drop helper column; keep karyotype_raw and new flags
    select(-karyotype_clean)
}


library(dplyr)
library(stringr)

recode_eln_keep_all <- function(df) {
  df %>%
    mutate(
      # start from the current eln column (whatever its name/case is)
      eln_raw = trimws(as.character(eln)),
      
      eln_clean = case_when(
        eln_raw == "Favorable"               ~ "Favorable",
        eln_raw == "Intermediate"            ~ "Intermediate",
        eln_raw == "Adverse"                 ~ "Adverse",
        eln_raw == "FavorableOrIntermediate" ~ "Intermediate",
        eln_raw == "IntermediateOrAdverse"   ~ "Intermediate",
        eln_raw %in% c("NonInitial", "NonAML") ~ "Unknown",
        TRUE                                 ~ "Unknown"   # catch anything odd
      ),
      
      eln_clean = factor(
        eln_clean,
        levels = c("Favorable", "Intermediate", "Adverse", "Unknown"),
        ordered = FALSE
      )
    )
}



## ============================================================
## 1) Apply to training set (c1)
## ============================================================

train <- read.csv("aml_baseline_clinical_beats_BM_c1.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

train2 <- add_cytogenetic_flags(train)

## Apply to training and testing
train2 <- recode_eln_keep_all(train2)

## Quick check
table(train2$eln_clean, useNA = "ifany")

train2 <- train2 %>%
  select(-karyotype_raw, -cyto_complex, -eln_raw, -eln)


write.csv(train2,
          "aml_baseline_clinical_beats_BM_c1_cyto.csv",
          row.names = FALSE)

## ============================================================
## 2) Apply to testing set (c2)
## ============================================================

test <- read.csv("aml_baseline_clinical_beats_BM_c2.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)

test2 <- add_cytogenetic_flags(test)

## Apply to training and testing
test2  <- recode_eln_keep_all(test2)

## Quick check
table(test2$eln_clean,  useNA = "ifany")

test2 <- test2 %>%
  select(-karyotype_raw, -cyto_complex, -eln_raw, -eln)


write.csv(test2,
          "aml_baseline_clinical_beats_BM_c2_cyto.csv",
          row.names = FALSE)
