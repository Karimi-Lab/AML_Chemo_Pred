###############################################################
## 0. Load libraries
###############################################################
library(dplyr)
library(stringr)
library(tibble)

###############################################################
## 1. Load the clinical RDA file
###############################################################

# Change path if needed
load("/Users/mehdi/Library/CloudStorage/OneDrive-King'sCollegeLondon(2)/Staffs/Sila_Gerlevic/AML/data/New_data/BEAT_BM_c2_clinical.Rda")

stopifnot(exists("clinical_df_c2_bm"))

df <- clinical_df_c2_bm


###############################################################
## 2. Helper recoding functions
###############################################################

# y/n fields → TRUE/FALSE
yn_to_logical <- function(x) {
  case_when(
    is.na(x) ~ NA,
    x %in% c("y", "Y", "yes", "Yes", "TRUE", "True") ~ TRUE,
    x %in% c("n", "N", "no", "No", "FALSE", "False") ~ FALSE,
    TRUE ~ NA
  )
}

# mutation calls like "positive"/"negative"
mut_to_logical <- function(x) {
  case_when(
    is.na(x) ~ NA,
    str_detect(tolower(x), "pos") ~ TRUE,
    str_detect(tolower(x), "neg") ~ FALSE,
    TRUE ~ NA
  )
}

# integer-coded mutations (0/1)
mut01 <- function(x) {
  case_when(
    is.na(x) ~ NA,
    x == 1 ~ TRUE,
    x == 0 ~ FALSE,
    TRUE ~ NA
  )
}


###############################################################
## 3. Clean mutation columns (0/1 fields)
###############################################################

mut01_cols <- c(
  "SRSF2","BRAF","WT1","CEBPA","IKZF1","NOTCH1","FLT3","FLT3-D835","JAK2","KIT",
  "DNMT3A","IDH1","PTPN11","CREBBP","SF3B1","BCOR","CCND2","GATA2","MLL2","SETBP1",
  "EZH2","TET2","ETV6","FBXW7","KRAS","NRAS","IDH2","CSF3R","KDM6A","CCND3",
  "CHEK2","CBL","U2AF1","MLL","SUZ12","FLT3-N676K","EP300","FAM5C","NOTCH2",
  "PHF6","JAK1","GATA1","ATM","MUTYH","STAG2","CIITA","CUX1","ZRSR2","FLT3-V491L",
  "POT1","MPL","CTCF","SMC1A","JAK3","BRCA2","NF1","BCORL1","PAX5","TCL1A","SOCS1",
  "RAD21","ASXL2","PDGFRB","CALR","BCL6","DNMT3A(p.R882H; MAF 16%)",
  "MYH11","CRLF2","KMT2D","FLT3-TKD","GNA13","CBLB","BCL2","CARD11",
  "No clinically significant variants identified","FAT4","RBBP6","ATG2B",
  "MED12","RTEL1","SH2B3","STAT3","CCND1","PRDM1","SYNE1","USH2A","CHD2",
  "PLCG2","DOCK8","FAT1","BLM","MSH2","PCLO","PRPF40B","TCF3"
)

mut01_cols <- mut01_cols[mut01_cols %in% colnames(df)]


###############################################################
## 4. Build AML baseline dataframe
###############################################################

aml_baseline <- df %>%
  mutate(
    ### Demographics
    sex = consensus_sex,
    age_dx = ageAtDiagnosis,
    
    ### Disease history
    prior_mds = yn_to_logical(priorMDS),
    prior_mpn = yn_to_logical(priorMPN),
    secondary_aml = prior_mds | prior_mpn,
    
    ### Response variables
    response_induction = responseToInductionTx,
    type_induction = typeInductionTx,
    
    ### Labs
    wbc = wbcCount,
    blasts_pb = as.numeric(`%.Blasts.in.PB`),
    blasts_bm = as.numeric(`%.Blasts.in.BM`),
    mono_pct = as.numeric(`%.Monocytes.in.PB`),
    lymph_pct = as.numeric(`%.Lymphocytes.in.PB`),
    platelet = plateletCount,
    hemoglobin = hemoglobin,
    hematocrit = hematocrit,
    albumin = albumin,
    creatinine = creatinine,
    ldh = LDH,
    
    ### Cytogenetics
    eln = ELN2017,
    karyotype = karyotype,
    
    ### FLT3-ITD (text)
    flt3_itd = mut_to_logical(`FLT3-ITD`),
    flt3_itd_ar = allelic_ratio,
    
    ### Key mutations (text or NA)
    npm1 = mut_to_logical(NPM1),
    runx1_text = RUNX1,   # needs parsing; complex strings
    asxl1_text = ASXL1,
    tp53_text = TP53
  ) %>%
  ### Parse RUNX1 / ASXL1 / TP53 complex entries → TRUE/FALSE
  mutate(
    runx1 = ifelse(!is.na(runx1_text) & runx1_text != "" & runx1_text != "NA", TRUE, FALSE),
    asxl1 = ifelse(!is.na(asxl1_text) & asxl1_text != "" & asxl1_text != "NA", TRUE, FALSE),
    tp53 = ifelse(!is.na(tp53_text) & tp53_text != "" & tp53_text != "NA", TRUE, FALSE)
  ) %>%
  
  ### Mutation columns coded as 0/1
  mutate(across(all_of(mut01_cols), mut01)) %>%
  
  ### Composite signatures
  mutate(
    spliceosome = SRSF2 | SF3B1 | U2AF1,
    cohesin = STAG2 | RAD21 | SMC1A,
    ras_pathway = KRAS | NRAS | PTPN11
  ) %>%
  
  ### Select the columns we want
  select(
    dbgap_rnaseq_sample, dbgap_subject_id,
    sex, age_dx, prior_mds, prior_mpn, secondary_aml,
    response_induction, type_induction,
    wbc, blasts_pb, blasts_bm, mono_pct, lymph_pct,
    platelet, hemoglobin, hematocrit, albumin, creatinine, ldh,
    eln, karyotype,
    flt3_itd, flt3_itd_ar,
    npm1, runx1, asxl1, tp53,
    DNMT3A, IDH1, IDH2, KRAS, NRAS, PTPN11,
    spliceosome, cohesin, ras_pathway
  )


###############################################################
## 5. Save results
###############################################################

aml_baseline$response_induction <- gsub(" i","", aml_baseline$response_induction)

library(dplyr)

aml_baseline_c2 <- aml_baseline %>%
  filter(
    response_induction %in% c("Refractory", "Complete Response")
  )

saveRDS(aml_baseline_c2, "aml_baseline_clinical_beats_BM_c2.rds")
write.csv(aml_baseline_c2, "aml_baseline_clinical_beats_BM_c2.csv", row.names = FALSE)



###############################################################
##  Done
###############################################################

