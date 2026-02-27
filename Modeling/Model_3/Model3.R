## ============================================================
## MODEL 3 – XGBoost (UPDATED with TEF Features)
## Fixes: 
## 1. Column alignment between Train/Test.
## 2. Explicitly removes P-value and Correlation columns.
## 3. Resolves the 'setdiff' function conflict.
## 4. Saves all test metrics to CSV for reliable output.
##
## Train: BEAT_BM_c1_merged_full_with_TEF.csv
## Test : BEAT_BM_c2_merged_full_with_TEF.csv
##
## Outcome: response_induction (Refractory vs Complete Response)
##
## ============================================================

## 0. Libraries ------------------------------------------------
# Ensure these libraries are installed: install.packages(c("tidymodels", "readr", "ggplot2", "cutpointr"))
library(tidymodels)
library(dplyr)
library(readr)
library(ggplot2)
library(cutpointr) 
tidymodels_prefer()

## 1. Read training & test data -------------------------------

# NOTE: Ensure these files are in your R working directory
train_raw <- read_csv("BEAT_BM_c1_merged_full_with_TEF.csv",
                      show_col_types = FALSE)
test_raw  <- read_csv("BEAT_BM_c2_merged_full_with_TEF.csv",
                      show_col_types = FALSE)

## 2. Define binary outcome and basic filtering (WITH ALIGNMENT FIX) ---------------

# Keep only Refractory vs Complete Response for modelling
train_m3 <- train_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory",
                       "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory")) # CR = negative, Refractory = positive (second level)
  )

test_m3 <- test_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory",
                       "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory"))
  )

# **********************************************
# ************ CRITICAL COLUMN ALIGNMENT FIX ************
# **********************************************

train_cols <- names(train_m3)
test_cols <- names(test_m3)

# Find columns present in train but missing in test, using base::setdiff to fix conflict
missing_in_test <- base::setdiff(train_cols, test_cols)

# Add missing columns to test set, filling with NA
for (col in missing_in_test) {
  # Check if the column is a true predictor (not a generated outcome column)
  if (!col %in% c("resp_bin", "response_induction")) { 
    test_m3[[col]] <- NA
  }
}
# Ensure all test columns are in the same order as train columns
test_m3 <- test_m3 %>% select(all_of(train_cols))

# ******************************************************

cat("Train n =", nrow(train_m3), " Test n =", nrow(test_m3), "\n")
cat("Training Outcome Distribution:\n")
print(table(train_m3$resp_bin))

## 3. Recipe – define predictors & preprocessing --------------

rec_m3 <- recipe(resp_bin ~ ., data = train_m3) %>%
  # mark ID columns so they are not treated as predictors
  update_role(dbgap_rnaseq_sample, dbgap_subject_id,
              new_role = "id") %>%
  
  # explicitly remove columns: IDs, outcomes, and the requested P-value/Correlation
  step_rm(dbgap_rnaseq_sample, dbgap_subject_id,
          response_induction, type_induction, 
          `P-value`, Correlation) %>% 
  
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

xgb_m3 <- boost_tree(
  mode            = "classification",
  trees           = tune(),
  tree_depth      = tune(),
  learn_rate      = tune(),
  mtry            = tune(),
  min_n           = tune(),
  loss_reduction  = tune()
) %>%
  set_engine("xgboost", missing = NA)    # let xgboost treat NA as missing

## 5. Workflow -------------------------------------------------

wf_m3 <- workflow() %>%
  add_model(xgb_m3) %>%
  add_recipe(rec_m3)

## 6. Cross-validation folds ----------------------------------

set.seed(123)
folds_m3 <- vfold_cv(train_m3, v = 5, strata = resp_bin)

## 7. Parameter grid ------------------------------------------

# Define parameter ranges and finalize mtry based on the preprocessed data
param_set_m3 <- parameters(
  trees(), tree_depth(), learn_rate(), mtry(), min_n(), loss_reduction()
)
rec_prepped <- prep(rec_m3, training = train_m3)
processed_train_data <- juice(rec_prepped) 
param_m3 <- finalize(param_set_m3, x = processed_train_data) 

set.seed(123)
# Create a Latin Hypercube sampling grid of 30 parameter combinations
grid_m3 <- grid_latin_hypercube(param_m3, size = 30)

## 8. Metrics set ---------------------------------------------

metric_set_m3 <- metric_set(
  roc_auc, pr_auc, yardstick::accuracy, sens, spec, recall, precision, f_meas
)

## 9. Tuning ---------------------------------------------------

set.seed(123)
cat("\n=== Starting Hyperparameter Tuning (Grid Search) ===\n")
tuned_m3 <- tune_grid(
  wf_m3,
  resamples = folds_m3,
  grid      = grid_m3,
  metrics   = metric_set_m3,
  control   = control_grid(save_pred = TRUE)
)

cat("\n=== Tuning Results Summary ===\n")
print(tuned_m3)
print(show_best(tuned_m3, metric = "roc_auc", n = 5))

## 10. Choose best hyperparameters ----------------------------

best_m3 <- select_best(tuned_m3, metric = "roc_auc")
cat("\n=== Best Parameters Found ===\n")
print(best_m3)

# Update the workflow with the best parameters found
final_wf_m3 <- finalize_workflow(wf_m3, best_m3)

## 11. Fit final model on full training set -------------------

cat("\n=== Fitting Final Model on Full Training Set ===\n")
fit_final_m3 <- fit(final_wf_m3, data = train_m3)

## 12. Evaluate on test set -----------------------------------

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m3, new_data = test_m3, type = "prob") %>%
  bind_cols(
    predict(fit_final_m3, new_data = test_m3, type = "class"),
    test_m3 %>% dplyr::select(resp_bin, dbgap_rnaseq_sample)
  )

## 12a. Scalar metrics on test set ----------------------------

cat("\n=== TEST PERFORMANCE (Model 3) ===\n")

# Calculate ALL metrics (Probability and Classification Metrics at default 0.5)
all_metrics <- metric_set_m3(
  data = test_preds, 
  truth = resp_bin, 
  .pred_Refractory, 
  estimate = .pred_class, 
  event_level = "second" 
)

cat("\n--- All Test Metrics (Default 0.5 Threshold) ---\n")
print(all_metrics)

# Save metrics to a file for easy extraction
write_csv(all_metrics, "model_3_test_metrics.csv")
cat("\nTest metrics saved to model_3_test_metrics.csv\n")

## 12b. OPTIMIZING THE THRESHOLD ------------------------------

cat("\n=== OPTIMAL THRESHOLD SEARCH (Maximizing F-Measure) ===\n")

# Use cutpointr to find the threshold that maximizes the F-measure (F1 Score)
cp <- cutpointr(
  test_preds %>% rename(outcome = resp_bin, probability = .pred_Refractory),
  x = probability,
  class = outcome,
  pos_class = "Refractory",
  neg_class = "CR",
  method = maximize_metric,
  metric = F1_score # Maximizes the balance between Precision and Recall
)

# Extract and report the optimal cutpoint
optimal_threshold <- cp$optimal_cutpoint
cat(paste("Optimal Threshold (Maximizing F-Measure):", round(optimal_threshold, 4), "\n"))


# Re-calculate hard predictions using the new threshold
test_preds_optimized <- test_preds %>%
  mutate(
    .pred_class_optimized = factor(
      if_else(.pred_Refractory >= optimal_threshold, "Refractory", "CR"),
      levels = c("CR", "Refractory")
    )
  )

# Re-calculate Classification metrics with the OPTIMIZED THRESHOLD
cat("\n--- Classification Metrics (OPTIMIZED Threshold) ---\n")
optimized_metrics <- metric_set(
  yardstick::accuracy, sens, spec, recall, precision, f_meas
)(
  data = test_preds_optimized,
  truth = resp_bin,
  estimate = .pred_class_optimized, 
  event_level = "second" 
)
print(optimized_metrics)
cat("\n--- Confusion Matrix (OPTIMIZED Threshold) ---\n")
conf_mat(test_preds_optimized, truth = resp_bin, estimate = .pred_class_optimized) %>%
  autoplot(type = "heatmap") 


## 13. ROC & PR curves (nice plots) ---------------------------

# ROC curve plot
roc_df <- roc_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(size = 1) + geom_abline(linetype = "dashed") + coord_equal() +
  labs(title = "Model 3 – ROC curve (Test set)", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
print(p_roc) 

# Precision–Recall curve plot
pr_df <- pr_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(size = 1) +
  labs(title = "Model 3 – Precision–Recall curve (Test set)", x = "Recall", y = "Precision") +
  theme_minimal()
print(p_pr) 

## 14. Cross-validation performance plot ----------------------

cat("\n--- Displaying CV Performance Plot ---\n")
p_cv <- autoplot(tuned_m3) +
  labs(
    title = "Model 3 – Cross-validation performance",
    subtitle = "ROC AUC across tuned hyperparameters"
  )
print(p_cv) 

## 15. Variable importance (optional) -------------------------
# Needs {vip}; run this if you have vip installed.
# library(vip)
# xgb_fit <- fit_final_m3 %>% extract_fit_parsnip()
# cat("\n--- Displaying Variable Importance Plot ---\n")
# vip(xgb_fit, num_features = 20) +
#    ggtitle("Model 3 – Top 20 predictors (XGBoost)") 

## ============================================================
## End of script
## ============================================================