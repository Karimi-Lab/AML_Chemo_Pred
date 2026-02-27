## ============================================================
## MODEL 3 – XGBoost (UPDATED with TEF Features, SHAP, and DR)
## Fixes: 
## 1. Column alignment between Train/Test.
## 2. Explicitly removes P-value and Correlation columns.
## 3. Resolves the 'setdiff' function conflict.
## 4. Saves all test metrics to CSV for reliable output.
##
## ADDED SECTIONS (15 & 16):
## - SHAP Feature Importance Analysis (4 plots)
## - Dimensionality Reduction (PCA & UMAP)
##
## **CRITICAL FIXES MAINTAINED:**
## - Static SHAP plot (15b) uses empty circles (shape=1) and robust data handling.
## - p-value plotting in section 15d is robust.
## ============================================================

MODEL_ID <- 3
N_TOP_FEATURES <- 20 # Total number of features to highlight in SHAP/DR plots

## 0. Libraries ------------------------------------------------
# Ensure these libraries are installed:
library(tidymodels)
library(dplyr)
library(readr)
library(ggplot2)
library(cutpointr) 

# --- New Libraries for SHAP/DR ---
library(fastshap)       # For SHAP calculation
library(patchwork)      # For arranging plots
library(uwot)           # For UMAP
library(scales)         # For plot formatting
library(Matrix)         # For matrix operations (SHAP)
library(shapviz)        # For SHAP plots (colored version)
library(tidyr)          # For data reshaping
library(rstatix)        # For wilcox.test, p-value formatting, and position
library(ggsignif)       # For adding significance bars
library(ggbeeswarm)     # For static SHAP beeswarm plot
# ---------------------------------

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

## 11. Fit final model on full training set and prepare data for SHAP/DR -------------------

cat("\n=== Fitting Final Model on Full Training Set ===\n")
fit_final_m3 <- fit(final_wf_m3, data = train_m3)

# Prepare Processed Data for SHAP/DR (CRITICAL STEP)
rec_final_m3 <- fit_final_m3 %>% extract_recipe()
X_train_processed <- bake(rec_final_m3, new_data = train_m3, all_predictors())
X_test_processed <- bake(rec_final_m3, new_data = test_m3, all_predictors())

# Get combined outcome data for plotting
y_train <- train_m3 %>% dplyr::select(resp_bin) %>% mutate(Dataset = "Train")
y_test <- test_m3 %>% dplyr::select(resp_bin) %>% mutate(Dataset = "Test")
y_combined <- bind_rows(y_train, y_test)


## 12. Evaluate on test set -----------------------------------

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m3, new_data = test_m3, type = "prob") %>%
  bind_cols(
    predict(fit_final_m3, new_data = test_m3, type = "class"),
    test_m3 %>% dplyr::select(resp_bin, dbgap_rnaseq_sample)
  )

## 12a. Scalar metrics on test set ----------------------------

cat("\n=== TEST PERFORMANCE (Model", MODEL_ID, ") ===\n")

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

## ... Section 12. Evaluate on test set -----------------------------------

# Generate class and probability predictions on the test set
test_preds_m3 <- predict(fit_final_m3, new_data = test_m3, type = "prob") %>%
  bind_cols(
    predict(fit_final_m3, new_data = test_m3, type = "class"),
    test_m3 %>% dplyr::select(resp_bin, dbgap_rnaseq_sample)
  )

# *******************************************************
# ************ CRITICAL STEP: SAVE PREDICTIONS **********
# *******************************************************
test_preds_m3 %>%
  # Rename the probability column to clearly identify the model
  rename(prob_Refractory = .pred_Refractory) %>% 
  # Select only the true outcome, the probability, and a unique ID
  select(dbgap_rnaseq_sample, resp_bin, prob_Refractory) %>%
  # Save to a file specific to the model
  write_csv("predictions_Model3.csv")
# *******************************************************

## 12a. Scalar metrics on test set ----------------------------
# ... rest of the script continues

## 13. ROC & PR curves (nice plots) ---------------------------

# ROC curve plot
roc_df <- roc_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(size = 1) + geom_abline(linetype = "dashed") + coord_equal() +
  labs(title = paste("Model", MODEL_ID, "– ROC curve (Test set)"), x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
print(p_roc) 
# Precision–Recall curve plot
pr_df <- pr_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(size = 1) +
  labs(title = paste("Model", MODEL_ID, "– Precision–Recall curve (Test set)"), x = "Recall", y = "Precision") +
  theme_minimal()
print(p_pr) 
## 14. Cross-validation performance plot ----------------------

cat("\n--- Displaying CV Performance Plot ---\n")
p_cv <- autoplot(tuned_m3) +
  labs(
    title = paste("Model", MODEL_ID, "– Cross-validation performance"),
    subtitle = "ROC AUC across tuned hyperparameters"
  )
print(p_cv) 
## 15. SHAP Feature Importance Analysis (Four Plots) ------------------------

xgb_engine <- fit_final_m3 %>% extract_fit_engine()

# Prediction wrapper function for SHAP (returns positive class probability)
p_fun_processed <- function(object, newdata) {
  newdata_matrix <- as.matrix(newdata)
  raw_pred <- predict(object, newdata = newdata_matrix)
  prob_matrix <- t(matrix(raw_pred, nrow = 2)) 
  return(prob_matrix[, 2]) # Extract probability of the positive class ("Refractory")
}

X_train_processed_complete <- X_train_processed
y_train_complete <- y_train

N_SHAP_INPUT <- nrow(X_train_processed_complete)
cat(paste("\n--- INPUT SAMPLE COUNT FOR SHAP (WITH NAs): ", N_SHAP_INPUT, " samples ---\n"))

cat("--- Calculating SHAP values on Training Set (nsim=25) ---\n")
set.seed(999) 
# Pass the feature data as a matrix with NAs.
shap_train <- fastshap::explain(
  xgb_engine, X = as.matrix(X_train_processed_complete), pred_wrapper = p_fun_processed, nsim = 25
)

N_SHAP_OUTPUT <- nrow(shap_train)
cat(paste("--- SHAP OUTPUT SAMPLE COUNT: ", N_SHAP_OUTPUT, " samples ---\n"))

# Guaranteed Fix: Match rows in case of silent failure
if (N_SHAP_OUTPUT < N_SHAP_INPUT) {
  warning(paste("SHAP output has fewer rows than input (", N_SHAP_OUTPUT, " vs ", N_SHAP_INPUT, 
                "). Subsetting X and Y to match SHAP output for plotting."))
  X_train_processed_complete <- X_train_processed_complete[1:N_SHAP_OUTPUT, ]
  y_train_complete <- y_train_complete[1:N_SHAP_OUTPUT, ]
}

# Create the shapviz object for the standard colored plot (15c)
shv <- shapviz(shap_train, X = X_train_processed_complete)

# Calculate Mean Absolute SHAP values and determine top features
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

top_features <- mean_abs_shap_train$Feature 

cat("\n--- Top", N_TOP_FEATURES, "Features by Mean Absolute SHAP ---\n")
print(mean_abs_shap_train)

# --- 15a. SHAP BAR CHART (Mean Absolute SHAP) ---
p_shap_barchart <- ggplot(mean_abs_shap_train, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(
    title = paste("Model", MODEL_ID, "– Global Feature Importance (Top", N_TOP_FEATURES, ")"),
    subtitle = "Ranked by Mean Absolute SHAP Value",
    x = "Mean Absolute SHAP Value",
    y = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

cat("\n--- Displaying SHAP Bar Chart ---\n")
print(p_shap_barchart) 

# Extract SHAP data for static ggplot beeswarm (avoids unexported shapviz::sv_data)
shap_data_for_plot <- shap_train %>%
  as.data.frame() %>%
  mutate(obs_id = row_number()) %>%
  tidyr::pivot_longer(
    cols = -obs_id,
    names_to = "feature",
    values_to = "value"
  ) %>%
  filter(feature %in% top_features)

# --- 15b. SHAP Signed Summary Plot (Static Color, EMPTY CIRCLES) ---
p_shap_signed_static <- ggplot(shap_data_for_plot, aes(x = value, y = feature)) +
  # Add a vertical line at SHAP = 0 (Neutral impact)
  geom_vline(xintercept = 0, color = "gray80", linetype = "dashed") + 
  # Use geom_quasirandom for the beeswarm effect.
  # Using shape = 1 for empty circle, color = "black" for the outline.
  ggbeeswarm::geom_quasirandom(shape = 1, color = "black", size = 2.5, alpha = 1, method = "pseudorandom") + 
  scale_y_discrete(limits = rev(top_features)) +
  labs(
    title = paste("Model", MODEL_ID, "– SHAP Signed Summary (Empty Circles)"),
    subtitle = "Dot position shows feature impact on prediction (right = Refractory).",
    x = "SHAP Value",
    y = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

cat("\n--- Displaying SHAP Signed Summary (Empty Circles) Plot ---\n")
print(p_shap_signed_static) 

# --- 15c. SHAP Signed Summary Plot (Colored by Feature Value) ---
p_shap_signed_colored <- sv_importance(shv, kind = "beeswarm", num_features = N_TOP_FEATURES) +
  labs(
    title = paste("Model", MODEL_ID, "– SHAP Signed Feature Summary (Colored by Value)"),
    subtitle = "Color indicates actual feature value. High values on the right push towards Refractory."
  ) +
  theme_minimal()

cat("\n--- Displaying SHAP Signed Summary (Colored by Value) Plot ---\n")
print(p_shap_signed_colored) 

## 15d. SHAP Distribution by Outcome (CR vs Refractory) with P-values ------------------------

cat("\n--- Performing Wilcoxon Rank-Sum Tests on Top SHAP Values ---\n")

# 1. Combine SHAP values with the true outcome (y_train_complete)
shap_long <- shap_train %>% 
  as.data.frame() %>%
  select(all_of(top_features)) %>%
  mutate(Outcome = y_train_complete$resp_bin) %>%
  tidyr::pivot_longer(
    cols = -Outcome,
    names_to = "Feature",
    values_to = "SHAP_Value"
  ) %>%
  # Ensure order matches importance plot
  mutate(Feature = factor(Feature, levels = rev(top_features))) 

# 2. Calculate P-values (Wilcoxon Rank-Sum Test)
p_values_df <- shap_long %>%
  group_by(Feature) %>%
  rstatix::wilcox_test(SHAP_Value ~ Outcome) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  rstatix::add_significance("p.adj") %>% 
  rstatix::add_xy_position(x = "Outcome", dodge = 0.8) %>% 
  mutate(
    p_text = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ paste("p.adj=", format.pval(p.adj, digits = 2))
    )
  )

cat("\n--- Top", N_TOP_FEATURES, "Features: SHAP P-values (CR vs Refractory) ---\n")
print(p_values_df %>% dplyr::select(Feature, p, p.adj, p.adj.signif))

# 3. Create the Visualization (Individual Paired Box-and-Jitter Plots)
cat("\n--- Displaying Paired SHAP Distribution Plots ---\n")
# Loop through the top features and generate a separate plot for each
plot_list <- list()
for (feature_name in top_features) {
  
  current_data <- shap_long %>% filter(Feature == feature_name)
  current_p_value <- p_values_df %>% filter(Feature == feature_name)
  
  max_y_data <- max(current_data$SHAP_Value)
  range_y <- max_y_data - min(current_data$SHAP_Value)
  sig_y_position <- max_y_data + (range_y * 0.15) 
  
  current_p_value <- current_p_value %>% 
    mutate(y.position = sig_y_position)
  
  p <- ggplot(current_data, aes(x = Outcome, y = SHAP_Value, fill = Outcome)) +
    geom_violin(alpha = 0.6, width = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
    expand_limits(y = sig_y_position * 1.05) + 
    geom_signif(
      data = current_p_value,
      aes(xmin = xmin, xmax = xmax, annotations = p_text, y = y.position),
      manual = TRUE,
      map_signif_level = FALSE,
      tip_length = 0.01,
      size = 0.5,
      inherit.aes = FALSE 
    ) +
    labs(
      title = paste("Feature:", feature_name),
      y = "SHAP Value",
      x = "Outcome"
    ) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          axis.title.x = element_blank())
  
  plot_list[[feature_name]] <- p
}

# Arrange and print the plots (using patchwork)
p_final_grid <- wrap_plots(plot_list, ncol = 5) + 
  plot_annotation(
    title = paste("Model", MODEL_ID, "– Paired SHAP Distribution by Outcome (Top", N_TOP_FEATURES, "Features)"),
    subtitle = "Wilcoxon Rank-Sum Test (Bonferroni adjusted p-values shown for CR vs. Refractory)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
                  plot.subtitle = element_text(hjust = 0.5))
  )

print(p_final_grid) 

## 16. Dimensionality Reduction (PCA & UMAP) on Top Features (Targeted Imputation) ---

# 1. Prepare combined data using only the Top Features
X_combined_processed <- bind_rows(X_train_processed, X_test_processed)
X_top_features <- X_combined_processed %>% dplyr::select(all_of(top_features))

# Get combined outcome data (must match rows of X_top_features)
y_combined_dr <- y_combined

# --- Targeted Imputation for PCA/UMAP ---
cat("\n--- Performing Targeted Imputation for Dimensionality Reduction ---\n")

# Create a function for column-wise median imputation
impute_median_col <- function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  }
  return(x)
}

# Apply imputation to the combined Top Features matrix
X_top_features_imputed <- X_top_features %>% 
  mutate(across(everything(), impute_median_col))

cat(paste("Total samples used for PCA/UMAP (after imputation):", nrow(X_top_features_imputed), "\n"))


# --- PCA ---
cat("\n--- Performing PCA on Top Features ---\n")
pca_res <- prcomp(X_top_features_imputed, scale = TRUE)
pca_df <- data.frame(pca_res$x[, 1:2]) %>% bind_cols(y_combined_dr) 
pc_variance <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = resp_bin)) +
  geom_point(alpha = 0.7, size = 2) + facet_wrap(~ Dataset) +
  scale_color_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
  labs(
    title = paste("Model", MODEL_ID, "– PCA on Top", N_TOP_FEATURES, "Features (All Samples)"),
    subtitle = paste0("PC1 Explains ", scales::percent(pc_variance[1]), "; PC2 Explains ", scales::percent(pc_variance[2])),
    color = "Outcome"
  ) +
  theme_minimal() + theme(legend.position = "bottom")

cat("\n--- Displaying PCA Plot ---\n")
print(p_pca) 

# --- UMAP ---
cat("\n--- Performing UMAP on Top Features (May take a moment) ---\n")
set.seed(999) 
umap_res <- uwot::umap(X_top_features_imputed, n_components = 2, verbose = FALSE)
umap_df <- data.frame(UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2]) %>% bind_cols(y_combined_dr)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = resp_bin)) +
  geom_point(alpha = 0.7, size = 2) + facet_wrap(~ Dataset) +
  scale_color_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
  labs(
    title = paste("Model", MODEL_ID, "– UMAP on Top", N_TOP_FEATURES, "Features (All Samples)"),
    color = "Outcome"
  ) +
  theme_minimal() + theme(legend.position = "bottom")

cat("\n--- Displaying UMAP Plot ---\n")
print(p_umap) 



library(ggplot2)
library(scales)

# 1. SETUP COLORS AND THEME
nature_colors <- c("CR" = "#4477AA", "Refractory" = "#EE6677")

theme_nature_pca <- function() {
  theme_bw(base_size = 12) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.position = "bottom"
    )
}

# 2. PCA PLOT WITH Y-AXIS LIMITS
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = resp_bin)) +
  stat_ellipse(aes(fill = resp_bin), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  geom_point(aes(fill = resp_bin), alpha = 0.8, size = 2.5, stroke = 0.3, shape = 21, color = "white") +
  facet_wrap(~ Dataset) +
  scale_color_manual(values = nature_colors) +
  scale_fill_manual(values = nature_colors) +
  labs(
    x = paste0("PC1 (", scales::percent(pc_variance[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", scales::percent(pc_variance[2], accuracy = 0.1), ")"),
    color = "Outcome"
  ) +
  # --- APPLYING THE Y-LIMITS HERE ---
  coord_cartesian(ylim = c(-5, 5)) + 
  theme_nature_pca()

print(p_pca)





library(ggplot2)
library(ggbeeswarm)
library(dplyr)

# 1. PREPARE THE DATA (Renaming step)
# This replaces the internal technical names with publication-ready labels
shap_data_for_plot <- shap_data_for_plot %>%
  mutate(feature = case_when(
    feature == "age_dx" ~ "Age",
    feature == "eln_clean_Favorable" ~ "ELN_Favorable",
    feature == "npm1_TRUE." ~ "NPM1_Mutation",
    TRUE ~ feature  # Keep others as is
  ))

# Update the top_features vector to match the new names so the ordering remains correct
top_features_renamed <- case_when(
  top_features == "age_dx" ~ "Age",
  top_features == "eln_clean_Favorable" ~ "ELN_Favorable",
  top_features == "npm1_TRUE." ~ "NPM1_Mutation",
  TRUE ~ top_features
)

# 2. DEFINE NATURE-STYLE SHAP THEME
theme_nature_shap <- function() {
  theme_classic(base_size = 12) %+replace%
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey95", linewidth = 0.5),
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5))
    )
}

# 3. GENERATE THE PLOT
p_shap_signed_static <- ggplot(shap_data_for_plot, aes(x = value, y = feature)) +
  # Solid thin black line for zero-intercept
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  
  # Nature-standard Beeswarm: Empty circles (shape 1) with transparency
  ggbeeswarm::geom_quasirandom(
    shape = 1, 
    color = "black", 
    size = 1.8, 
    stroke = 0.4, 
    alpha = 0.7, 
    method = "pseudorandom"
  ) + 
  
  # Apply renamed and ordered features to Y-axis
  scale_y_discrete(limits = rev(top_features_renamed)) +
  
  labs(
    title = paste("Model 3 – SHAP Feature Impact"),
    x = "SHAP Value (Impact on Model Output)",
    y = NULL
  ) +
  
  theme_nature_shap()

# Display
print(p_shap_signed_static)

## 17. Variable importance (optional) -------------------------
# Needs {vip}; run this if you have vip installed.
# library(vip)
# xgb_fit <- fit_final_m3 %>% extract_fit_parsnip()
# cat("\n--- Displaying Variable Importance Plot ---\n")
# vip(xgb_fit, num_features = 20) +
#    ggtitle(paste("Model", MODEL_ID, "– Top 20 predictors (XGBoost)")) 

## ============================================================
## End of script
## ============================================================