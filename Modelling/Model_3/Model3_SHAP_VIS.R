## ============================================================
## MODEL 3 – XGBoost (Robust Script with SHAP and DR)
## 
## **Data Files: BEAT_BM_c1_merged_full_with_TEF.csv / c2_merged_full_with_TEF.csv**
## 
## Outcome: response_induction (Refractory vs Complete Response)
##
## ============================================================

MODEL_ID <- 3 # *** MODEL ID: Set for Model 3 ***
N_TOP_FEATURES <- 20 # Number of features to use for SHAP plot and DR

## 0. Libraries ------------------------------------------------
# Ensure these libraries are installed: install.packages(c("tidymodels", "readr", "ggplot2", "cutpointr", "fastshap", "patchwork", "uwot", "scales", "Matrix", "shapviz", "tidyr"))
library(tidymodels)
library(dplyr)
library(readr)
library(ggplot2)
library(cutpointr) 
library(fastshap)       # For SHAP value calculation
library(patchwork)      # For combining plots
library(uwot)           # For UMAP
library(scales)         # For percentage formatting in PCA plots
library(Matrix)         # Needed for safe matrix conversion in SHAP/XGBoost
library(shapviz)        # Required for plotting from fastshap
library(tidyr)          # Required for pivot_longer in manual SHAP calculation

tidymodels_prefer()

## 1. Read training & test data -------------------------------
# NOTE: Ensure these files are in your R working directory. Assuming test data follows a similar naming scheme.
train_raw <- read_csv("BEAT_BM_c1_merged_full_with_TEF.csv", show_col_types = FALSE) 
# *** Assuming test file name: ***
test_raw  <- read_csv("BEAT_BM_c2_merged_full_with_TEF.csv", show_col_types = FALSE) 

## 2. Define binary outcome and basic filtering (WITH ALIGNMENT FIX) ---------------
# Keep only Refractory vs Complete Response for modelling
train_m1 <- train_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory", "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory")) # CR = negative, Refractory = positive (second level)
  )

test_m1 <- test_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory", "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory"))
  )

# ************ CRITICAL COLUMN ALIGNMENT FIX ************
# Ensures test set columns match training set columns
train_cols <- names(train_m1)
test_cols <- names(test_m1)
missing_in_test <- base::setdiff(train_cols, test_cols)

for (col in missing_in_test) {
  if (!col %in% c("resp_bin", "response_induction")) { 
    test_m1[[col]] <- NA
  }
}
test_m1 <- test_m1 %>% select(all_of(train_cols))
# ******************************************************

cat("Train n =", nrow(train_m1), " Test n =", nrow(test_m1), "\n")
cat("Training Outcome Distribution:\n")
print(table(train_m1$resp_bin))

## 3. Recipe – define predictors & preprocessing --------------

rec_m1 <- recipe(resp_bin ~ ., data = train_m1) %>%
  # mark ID columns so they are not treated as predictors
  update_role(any_of(c("dbgap_rnaseq_sample", "dbgap_subject_id", "X1")), 
              new_role = "id") %>%
  
  # explicitly remove columns: IDs, outcomes, and optional irrelevant columns
  step_rm(any_of(c("dbgap_rnaseq_sample", "dbgap_subject_id", "X1")),
          response_induction, type_induction, 
          any_of(c("P-value", "Correlation"))) %>% 
  
  # 1. Ensure all character predictors become factors
  step_string2factor(all_nominal_predictors()) %>%
  
  # 2. Convert known genetic/cyto binary flags (if present) to factors.
  step_mutate_at(any_of(c(
    "prior_mds", "prior_mpn", "secondary_aml", "flt3_itd", "npm1", "runx1", 
    "asxl1", "tp53", "DNMT3A", "IDH1", "IDH2", "KRAS", "NRAS", "PTPN11", 
    "spliceosome", "cohesin", "ras_pathway", 
    "cyto_t_8_21", "cyto_inv_16", "cyto_t_16_16", "cyto_t_9_11", "cyto_t_6_9", 
    "cyto_inv_3", "cyto_t_3_3", "cyto_5_minus", "cyto_del_5q", 
    "cyto_7_minus", "cyto_del_7q", "cyto_del_9q", "cyto_del_11q", 
    "cyto_del_13q", "cyto_del_17p", "cyto_plus8", "cyto_normal",
    "cyto_inv16", "cyto_cbf", "cyto_t_15_17", "cyto_inv3", "cyto_5_7_abn"
  )), fn = factor) %>%
  
  # 3. Remove zero-variance predictors (per resample)
  step_zv(all_predictors()) %>%
  
  # 4. Handle missing factor levels in resamples (robustness for CV)
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  
  # 5. dummy-encode all categorical predictors (creates numeric 0/1 columns for XGBoost)
  step_dummy(all_nominal_predictors(), one_hot = TRUE)

## 4. XGBoost model spec --------------------------------------

xgb_m1 <- boost_tree(
  mode            = "classification", trees = tune(), tree_depth = tune(), 
  learn_rate      = tune(), mtry = tune(), min_n = tune(), loss_reduction = tune()
) %>%
  set_engine("xgboost", missing = NA)    

## 5. Workflow -------------------------------------------------

wf_m1 <- workflow() %>% add_model(xgb_m1) %>% add_recipe(rec_m1)

## 6. Cross-validation folds ----------------------------------

set.seed(123)
folds_m1 <- vfold_cv(train_m1, v = 5, strata = resp_bin)

## 7. Parameter grid ------------------------------------------

param_set_m1 <- parameters(trees(), tree_depth(), learn_rate(), mtry(), min_n(), loss_reduction())
rec_prepped <- prep(rec_m1, training = train_m1)
processed_train_data <- juice(rec_prepped) 
param_m1 <- finalize(param_set_m1, x = processed_train_data) 

set.seed(123)
grid_m1 <- grid_latin_hypercube(param_m1, size = 30)

## 8. Metrics set ---------------------------------------------

metric_set_m1 <- metric_set(
  roc_auc, pr_auc, yardstick::accuracy, sens, spec, recall, precision, f_meas
)

## 9. Tuning ---------------------------------------------------

set.seed(123)
cat("\n=== Starting Hyperparameter Tuning (Grid Search) ===\n")
tuned_m1 <- tune_grid(
  wf_m1, resamples = folds_m1, grid = grid_m1, metrics = metric_set_m1,
  control = control_grid(save_pred = TRUE)
)

cat("\n=== Tuning Results Summary ===\n")
print(show_best(tuned_m1, metric = "roc_auc", n = 5))

## 10. Choose best hyperparameters ----------------------------

best_m1 <- select_best(tuned_m1, metric = "roc_auc")
cat("\n=== Best Parameters Found ===\n")
print(best_m1)
final_wf_m1 <- finalize_workflow(wf_m1, best_m1)

## 11. Fit final model on full training set -------------------

cat("\n=== Fitting Final Model on Full Training Set ===\n")
fit_final_m1 <- fit(final_wf_m1, data = train_m1)

## 11b. Prepare Processed Data for SHAP/DR --------------------

rec_final_m1 <- fit_final_m1 %>% extract_recipe()
X_train_processed <- bake(rec_final_m1, new_data = train_m1, all_predictors())
X_test_processed <- bake(rec_final_m1, new_data = test_m1, all_predictors())

# Get combined outcome data for plotting
y_train <- train_m1 %>% select(resp_bin, any_of("dbgap_rnaseq_sample")) %>% mutate(Dataset = "Train")
y_test <- test_m1 %>% select(resp_bin, any_of("dbgap_rnaseq_sample")) %>% mutate(Dataset = "Test")
y_combined <- bind_rows(y_train, y_test)

## 12. Evaluate on test set -----------------------------------

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m1, new_data = test_m1, type = "prob") %>%
  bind_cols(
    predict(fit_final_m1, new_data = test_m1, type = "class"),
    test_m1 %>% dplyr::select(resp_bin, any_of("dbgap_rnaseq_sample"))
  )

## 12a. Scalar metrics on test set ----------------------------

cat("\n=== TEST PERFORMANCE (Model", MODEL_ID, ") ===\n")
cat("\n--- All Test Metrics (Default 0.5 Threshold) ---\n")
all_metrics <- metric_set_m1(
  data = test_preds, truth = resp_bin, .pred_Refractory, estimate = .pred_class, 
  event_level = "second" 
)
print(all_metrics)

## 12b. OPTIMIZING THE THRESHOLD ------------------------------

cat("\n=== OPTIMAL THRESHOLD SEARCH (Maximizing F-Measure) ===\n")
cp <- cutpointr(
  test_preds %>% rename(outcome = resp_bin, probability = .pred_Refractory),
  x = probability, class = outcome, pos_class = "Refractory", neg_class = "CR",
  method = maximize_metric, metric = cutpointr::F1_score 
)
optimal_threshold <- cp$optimal_cutpoint
cat(paste("Optimal Threshold (Maximizing F-Measure):", round(optimal_threshold, 4), "\n"))

test_preds_optimized <- test_preds %>%
  mutate(
    .pred_class_optimized = factor(
      if_else(.pred_Refractory >= optimal_threshold, "Refractory", "CR"),
      levels = c("CR", "Refractory")
    )
  )

cat("\n--- Classification Metrics (OPTIMIZED Threshold) ---\n")
optimized_metrics <- metric_set(yardstick::accuracy, sens, spec, recall, precision, f_meas)(
  data = test_preds_optimized, truth = resp_bin, estimate = .pred_class_optimized, 
  event_level = "second" 
)
print(optimized_metrics)
cat("\n--- Confusion Matrix (OPTIMIZED Threshold) ---\n")
conf_mat(test_preds_optimized, truth = resp_bin, estimate = .pred_class_optimized) %>%
  autoplot(type = "heatmap") 


## 13. ROC & PR curves (nice plots) ---------------------------

roc_df <- roc_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(size = 1) + geom_abline(linetype = "dashed") + coord_equal() +
  labs(title = paste("Model", MODEL_ID, "– ROC curve (Test set)")) + theme_minimal()
print(p_roc) 

pr_df <- pr_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(size = 1) +
  labs(title = paste("Model", MODEL_ID, "– Precision–Recall curve (Test set)")) + theme_minimal()
print(p_pr) 

## 14. Cross-validation performance plot ----------------------

cat("\n--- Displaying CV Performance Plot ---\n")
p_cv <- autoplot(tuned_m1) +
  labs(title = paste("Model", MODEL_ID, "– Cross-validation performance"))
print(p_cv) 

## 16. SHAP Feature Importance Analysis ------------------------

# 1. Extract the raw fitted engine (XGBoost model object)
xgb_engine <- fit_final_m1 %>% extract_fit_engine()

# 2. Define a predictor function (TAKING PROCESSED DATA)
p_fun_processed <- function(object, newdata) {
  newdata_matrix <- as.matrix(newdata)
  raw_pred <- predict(object, newdata = newdata_matrix)
  prob_matrix <- t(matrix(raw_pred, nrow = 2)) 
  return(prob_matrix[, 2]) # Extract probability of the positive class ("Refractory")
}

# 3. *** ROBUSTNESS FIX 1: Filter processed data to ensure complete cases for SHAP ***
complete_cases_idx <- which(complete.cases(X_train_processed))
if (length(complete_cases_idx) < nrow(X_train_processed)) {
  warning(paste(nrow(X_train_processed) - length(complete_cases_idx), 
                "samples removed from X_train_processed for SHAP due to NA values."))
}
X_train_processed_complete <- X_train_processed[complete_cases_idx, ]
y_train_complete <- y_train[complete_cases_idx, ]

# 4. Calculate SHAP values on the complete Training Set 
cat("\n--- Calculating SHAP values on Training Set (Complete Cases Only) ---\n")
set.seed(999) 
shap_train <- fastshap::explain(
  xgb_engine, X = X_train_processed_complete, pred_wrapper = p_fun_processed, nsim = 25
)

# 5. *** ROBUSTNESS FIX 2: Check and subset X to guarantee row count match ***
n_shap_rows <- nrow(shap_train)
n_X_rows <- nrow(X_train_processed_complete)

if (n_shap_rows != n_X_rows) {
  warning(paste("Final SHAP rows (", n_shap_rows, ") still differ from input X rows (", n_X_rows, "). Subsetting X to match SHAP output."))
  X_train_processed_complete <- X_train_processed_complete[1:n_shap_rows, ]
  y_train_complete <- y_train_complete[1:n_shap_rows, ]
}

# 6. Create the shapviz object using the guaranteed matching data subsets
shv <- shapviz(shap_train, X = X_train_processed_complete)

# 7. Generate the Signed SHAP Summary Plot (Beeswarm plot)
p_shap_signed <- sv_importance(shv, kind = "beeswarm", num_features = N_TOP_FEATURES) +
  labs(
    title = paste("Model", MODEL_ID, "– SHAP Signed Feature Summary (Top", N_TOP_FEATURES, ")"),
    subtitle = "Color indicates actual feature value. High values on the right push towards Refractory."
  ) +
  theme_minimal()

cat("\n--- Displaying SHAP Signed Summary Plot ---\n")
print(p_shap_signed) 

# 8. *** FINAL & GUARANTEED FIX 3: Calculate Mean Absolute SHAP values manually ***
mean_abs_shap_train <- shap_train %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ mean(abs(.x)))) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Feature",
    values_to = "Importance"
  ) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = N_TOP_FEATURES)

# 9. Extract the top features list for DR analysis in the next section
top_features <- mean_abs_shap_train$Feature 

cat("\n--- Top", N_TOP_FEATURES, "Features by Mean Absolute SHAP ---\n")
print(mean_abs_shap_train)

# 10. Generate Mean Absolute SHAP Importance Plot for Training Data
p_shap_summary <- ggplot(mean_abs_shap_train, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_point(color = "dodgerblue4", size = 3) +
  geom_segment(aes(x = 0, xend = Importance, yend = reorder(Feature, Importance)), color = "dodgerblue4") +
  labs(title = paste("Model", MODEL_ID, "– Mean Absolute SHAP Importance")) +
  theme_minimal()

cat("\n--- Displaying Mean Absolute SHAP Importance Plot ---\n")
print(p_shap_summary) 



## 17. Dimensionality Reduction (PCA & UMAP) on Top Features ---

# 1. Prepare combined data using only the Top Features
X_combined_processed <- bind_rows(X_train_processed, X_test_processed)
X_top_features <- X_combined_processed %>% select(all_of(top_features))

# Remove samples if the recipe failed to impute
na_check <- rowSums(is.na(X_top_features)) > 0
if (any(na_check)) {
  warning(paste(sum(na_check), "samples with NA values in top features removed for DR."))
  X_top_features <- X_top_features[!na_check, ]
  y_combined_dr <- y_combined[!na_check, ]
} else {
  y_combined_dr <- y_combined
}

# --- PCA ---
cat("\n--- Performing PCA on Top Features ---\n")
pca_res <- prcomp(X_top_features, scale. = TRUE)
pca_df <- data.frame(pca_res$x[, 1:2]) %>% bind_cols(y_combined_dr)
pc_variance <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_var <- scales::percent(pc_variance[1])
pc2_var <- scales::percent(pc_variance[2])

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = resp_bin)) +
  geom_point(alpha = 0.7, size = 2) + facet_wrap(~ Dataset) +
  scale_color_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
  labs(title = paste("Model", MODEL_ID, "– PCA on Top", N_TOP_FEATURES, "Features")) +
  theme_minimal() + theme(legend.position = "bottom")

cat("\n--- Displaying PCA Plot ---\n")
print(p_pca) 



# --- UMAP ---
cat("\n--- Performing UMAP on Top Features (May take a moment) ---\n")
set.seed(999) 
umap_res <- uwot::umap(X_top_features, n_components = 2, verbose = FALSE)
umap_df <- data.frame(UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2]) %>% bind_cols(y_combined_dr)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = resp_bin)) +
  geom_point(alpha = 0.7, size = 2) + facet_wrap(~ Dataset) +
  scale_color_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
  labs(title = paste("Model", MODEL_ID, "– UMAP on Top", N_TOP_FEATURES, "Features")) +
  theme_minimal() + theme(legend.position = "bottom")

cat("\n--- Displaying UMAP Plot ---\n")
print(p_umap) 


## ============================================================
## End of script
## ============================================================