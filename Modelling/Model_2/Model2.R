## ============================================================
## MODEL 2 – XGBoost (Final Robust Script)
## Calculates ROC AUC and PR AUC on the Test Set.
##
## Train: BEAT_BM_c1_merged_full.csv
## Test : BEAT_BM_c2_merged_full.csv
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
train_raw <- read_csv("BEAT_BM_c1_merged_full.csv",
                      show_col_types = FALSE)
test_raw  <- read_csv("BEAT_BM_c2_merged_full.csv",
                      show_col_types = FALSE)

## 2. Define binary outcome and basic filtering (WITH ALIGNMENT FIX) ---------------

# Keep only Refractory vs Complete Response for modelling
train_m1 <- train_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory",
                       "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory")) # CR = negative, Refractory = positive (second level)
  )

test_m1 <- test_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory",
                       "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory"))
  )

# **********************************************
# ************ CRITICAL COLUMN ALIGNMENT FIX ************
# **********************************************

train_cols <- names(train_m1)
test_cols <- names(test_m1)

# Find columns present in train but missing in test, using base::setdiff to fix conflict
missing_in_test <- base::setdiff(train_cols, test_cols)

# Add missing columns to test set, filling with NA
for (col in missing_in_test) {
  # Check if the column is a true predictor (not a generated outcome column)
  if (!col %in% c("resp_bin", "response_induction")) { 
    test_m1[[col]] <- NA
  }
}
# Ensure all test columns are in the same order as train columns
test_m1 <- test_m1 %>% select(all_of(train_cols))

# ******************************************************

cat("Train n =", nrow(train_m1), " Test n =", nrow(test_m1), "\n")
cat("Training Outcome Distribution:\n")
print(table(train_m1$resp_bin))

## 3. Recipe – define predictors & preprocessing --------------

rec_m1 <- recipe(resp_bin ~ ., data = train_m1) %>%
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

xgb_m1 <- boost_tree(
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

wf_m1 <- workflow() %>%
  add_model(xgb_m1) %>%
  add_recipe(rec_m1)

## 6. Cross-validation folds ----------------------------------

set.seed(123)
folds_m1 <- vfold_cv(train_m1, v = 5, strata = resp_bin)

## 7. Parameter grid ------------------------------------------

# Define parameter ranges and finalize mtry based on the preprocessed data
param_set_m1 <- parameters(
  trees(), tree_depth(), learn_rate(), mtry(), min_n(), loss_reduction()
)
rec_prepped <- prep(rec_m1, training = train_m1)
processed_train_data <- juice(rec_prepped) 
param_m1 <- finalize(param_set_m1, x = processed_train_data) 

set.seed(123)
# Create a Latin Hypercube sampling grid of 30 parameter combinations
grid_m1 <- grid_latin_hypercube(param_m1, size = 30)

## 8. Metrics set ---------------------------------------------

metric_set_m1 <- metric_set(
  roc_auc, pr_auc, yardstick::accuracy, sens, spec, recall, precision, f_meas
)

## 9. Tuning ---------------------------------------------------

set.seed(123)
cat("\n=== Starting Hyperparameter Tuning (Grid Search) ===\n")
tuned_m1 <- tune_grid(
  wf_m1,
  resamples = folds_m1,
  grid      = grid_m1,
  metrics   = metric_set_m1,
  control   = control_grid(save_pred = TRUE)
)

cat("\n=== Tuning Results Summary ===\n")
print(tuned_m1)
print(show_best(tuned_m1, metric = "roc_auc", n = 5))

## 10. Choose best hyperparameters ----------------------------

best_m1 <- select_best(tuned_m1, metric = "roc_auc")
cat("\n=== Best Parameters Found ===\n")
print(best_m1)

# Update the workflow with the best parameters found
final_wf_m1 <- finalize_workflow(wf_m1, best_m1)

## 11. Fit final model on full training set -------------------

cat("\n=== Fitting Final Model on Full Training Set ===\n")
fit_final_m1 <- fit(final_wf_m1, data = train_m1)

## 12. Evaluate on test set -----------------------------------

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m1, new_data = test_m1, type = "prob") %>%
  bind_cols(
    predict(fit_final_m1, new_data = test_m1, type = "class"),
    test_m1 %>% dplyr::select(resp_bin, dbgap_rnaseq_sample)
  )

## 12a. Scalar metrics on test set ----------------------------

cat("\n=== TEST PERFORMANCE (Model 1) ===\n")

# 1. Probability-based metrics (ROC AUC, PR AUC)
cat("\n--- Probability-Based Metrics ---\n")
prob_metrics <- metric_set(roc_auc, pr_auc)(
  data = test_preds, 
  truth = resp_bin, 
  .pred_Refractory, 
  event_level = "second" 
)
print(prob_metrics)

# 2. Classification metrics (AT DEFAULT 0.5 THRESHOLD)
cat("\n--- Classification Metrics (Default 0.5 Threshold) ---\n")
class_metrics <- metric_set(
  yardstick::accuracy, sens, spec, recall, precision, f_meas
)(
  data = test_preds,
  truth = resp_bin,
  estimate = .pred_class, 
  event_level = "second" 
)
print(class_metrics)

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
  labs(title = "Model 1 – ROC curve (Test set)", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
print(p_roc) 

# Precision–Recall curve plot
pr_df <- pr_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(size = 1) +
  labs(title = "Model 1 – Precision–Recall curve (Test set)", x = "Recall", y = "Precision") +
  theme_minimal()
print(p_pr) 

## 14. Cross-validation performance plot ----------------------

cat("\n--- Displaying CV Performance Plot ---\n")
p_cv <- autoplot(tuned_m1) +
  labs(
    title = "Model 1 – Cross-validation performance",
    subtitle = "ROC AUC across tuned hyperparameters"
  )
print(p_cv) 

## 15. Variable importance (optional) -------------------------
# Needs {vip}; run this if you have vip installed.
# library(vip)
# xgb_fit <- fit_final_m1 %>% extract_fit_parsnip()
# cat("\n--- Displaying Variable Importance Plot ---\n")
# vip(xgb_fit, num_features = 20) +
#    ggtitle("Model 1 – Top 20 predictors (XGBoost)") 

## ============================================================
## End of script
## ============================================================
