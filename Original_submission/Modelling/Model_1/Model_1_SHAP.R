## ============================================================
## MODEL 1 – FINAL INTEGRATED SCRIPT (XGBoost + SHAP + DR)
## 
## **Final Integrated Interpretability Views (4 SHAP Plots):**
## 1. Global Feature Importance (SHAP Bar Chart)
## 2. SHAP Signed Summary (Static Color - Empty Circles) 
## 3. SHAP Signed Summary (Colored by Feature Value - PREVIOUS)
## 4. Paired SHAP Distribution (Violin/Box Plots) by Outcome with visible adjusted P-values.
## 
## **DEFINITIVE ROBUST FIXES APPLIED:**
## - CRITICAL: Removed reliance on unexported shapviz::sv_data() by manually converting SHAP matrix to long format (Section 9b).
## - FIXED: SHAP Signed Summary (9b) now uses empty, plain black circles (shape=1).
## - Corrected rstatix calls, geom_signif inheritance, and P-value axis expansion (Section 9d).
## - Corrected typo 'data.data.frame' to 'data.frame' (Section 10).
## 
## ============================================================

MODEL_ID <- 1 
N_TOP_FEATURES <- 20 # Total number of features to highlight in SHAP/DR plots
DEP_PLOT_FEATURE <- "age" # Example feature for dependency plot 

## 0. Libraries ------------------------------------------------
# Ensure these libraries are installed: tidymodels, readr, ggplot2, cutpointr, fastshap, patchwork, uwot, scales, Matrix, shapviz, tidyr, rstatix, ggsignif, ggbeeswarm
library(tidymodels)
library(dplyr)
library(readr)
library(ggplot2)
library(cutpointr) 
library(fastshap)       
library(patchwork)      
library(uwot)           
library(scales)         
library(Matrix)         
library(shapviz)        
library(tidyr)        
library(rstatix)      # For wilcox.test, p-value formatting, and position
library(ggsignif)     # For adding significance bars to plots
library(ggbeeswarm)   # For static SHAP beeswarm plot

tidymodels_prefer()

## 1. Read training & test data -------------------------------

# NOTE: Ensure these files are in your R working directory
train_raw <- read_csv("aml_baseline_clinical_beats_BM_c1_cyto.csv", show_col_types = FALSE)
test_raw <- read_csv("aml_baseline_clinical_beats_BM_c2_cyto.csv", show_col_types = FALSE)

# 1a. Define binary outcome and filter
train_m1 <- train_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory", "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory")) 
  )

test_m1 <- test_raw %>%
  filter(response_induction %in% c("Refractory", "Complete Response")) %>%
  mutate(
    resp_bin = if_else(response_induction == "Refractory", "Refractory", "CR"),
    resp_bin = factor(resp_bin, levels = c("CR", "Refractory")) 
  )

cat("Train n =", nrow(train_m1), " Test n =", nrow(test_m1), "\n")
print(table(train_m1$resp_bin))

## 2. Recipe – define predictors & preprocessing (NA Preserved for XGBoost) --------------

rec_m1 <- recipe(resp_bin ~ ., data = train_m1) %>%
  # mark ID columns
  update_role(any_of(c("dbgap_subject_id", "X1", "dbgap_rnaseq_sample")), new_role = "id") %>%
  
  # remove irrelevant columns
  step_rm(any_of(c("dbgap_subject_id", "X1", "dbgap_rnaseq_sample")),
          response_induction, type_induction) %>% 
  
  # 1. Ensure all character predictors become factors
  step_string2factor(all_nominal_predictors()) %>%
  
  # 2. CRITICAL FIX: Explicitly convert all known gene/disease flags to factors.
  step_mutate_at(any_of(c(
    "gender", 
    "prior_mds", "prior_mpn", "secondary_aml", 
    "flt3_itd", "npm1", "runx1", "asxl1", "tp53", "DNMT3A", "IDH1", "IDH2", 
    "KRAS", "NRAS", "PTPN11", "spliceosome", "cohesin", "ras_pathway", 
    "cyto_t_8_21", "cyto_inv_16", "cyto_t_16_16", "cyto_t_9_11", "cyto_t_6_9", 
    "cyto_inv_3", "cyto_t_3_3", "cyto_5_minus", "cyto_del_5q", 
    "cyto_7_minus", "cyto_del_7q", "cyto_del_9q", "cyto_del_11q", 
    "cyto_del_13q", "cyto_del_17p", "cyto_plus8", "cyto_normal",
    "cyto_inv16", "cyto_cbf", "cyto_t_15_17", "cyto_inv3", "cyto_5_7_abn"
  )), fn = factor) %>%
  
  # 3. Preserve Numeric NAs (no imputation for modeling)
  
  # 4. Remove zero-variance predictors
  step_zv(all_predictors()) %>%
  
  # 5. Handle missing factor levels (for robustness)
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors(), new_level = "missing_category") %>%
  
  # 6. dummy-encode all categorical predictors 
  step_dummy(all_nominal_predictors(), one_hot = TRUE)

## 3. XGBoost model spec --------------------------------------

xgb_m1 <- boost_tree(
  mode            = "classification", trees = tune(), tree_depth = tune(), 
  learn_rate      = tune(), mtry = tune(), min_n = tune(), loss_reduction = tune()
) %>%
  set_engine("xgboost", missing = NA) 

## 4. Workflow, Folds, and Grid Setup --------------------------

wf_m1 <- workflow() %>% add_model(xgb_m1) %>% add_recipe(rec_m1)

set.seed(123)
folds_m1 <- vfold_cv(train_m1, v = 5, strata = resp_bin)

param_set_m1 <- parameters(trees(), tree_depth(), learn_rate(), mtry(), min_n(), loss_reduction())
rec_prepped <- prep(rec_m1, training = train_m1)
processed_train_data <- juice(rec_prepped) 
param_m1 <- finalize(param_set_m1, x = processed_train_data) 

set.seed(123)
grid_m1 <- grid_latin_hypercube(param_m1, size = 30)

## 5. Metrics set and Tuning ----------------------------------

metric_set_m1 <- metric_set(
  roc_auc, pr_auc, yardstick::accuracy, sens, spec, recall, precision, f_meas
)

set.seed(123)
cat("\n=== Starting Hyperparameter Tuning (Grid Search) ===\n")
tuned_m1 <- tune_grid(
  wf_m1, resamples = folds_m1, grid = grid_m1, metrics = metric_set_m1,
  control = control_grid(save_pred = TRUE)
)

cat("\n=== Tuning Results Summary ===\n")
print(show_best(tuned_m1, metric = "roc_auc", n = 5))

## 6. Final Fit and Processed Data Preparation -------------------

best_m1 <- select_best(tuned_m1, metric = "roc_auc")
cat("\n=== Best Parameters Found ===\n")
print(best_m1)
final_wf_m1 <- finalize_workflow(wf_m1, best_m1)

cat("\n=== Fitting Final Model on Full Training Set ===\n")
fit_final_m1 <- fit(final_wf_m1, data = train_m1)

# Prepare Processed Data for SHAP/DR
rec_final_m1 <- fit_final_m1 %>% extract_recipe()
X_train_processed <- bake(rec_final_m1, new_data = train_m1, all_predictors())
X_test_processed <- bake(rec_final_m1, new_data = test_m1, all_predictors())

# Get combined outcome data for plotting
y_train <- train_m1 %>% select(resp_bin) %>% mutate(Dataset = "Train")
y_test <- test_m1 %>% select(resp_bin) %>% mutate(Dataset = "Test")
y_combined <- bind_rows(y_train, y_test)

## 7. Evaluate on test set (Sections 7a, 7b, and 8 included here) -----------------------------------

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m1, new_data = test_m1, type = "prob") %>%
  bind_cols(
    predict(fit_final_m1, new_data = test_m1, type = "class"),
    test_m1 %>% dplyr::select(resp_bin, dbgap_subject_id) # <-- Make sure dbgap_subject_id is here!
  )


## 7. Evaluate on test set (Sections 7a, 7b, and 8 included here) -----------------------------------

# # Generate class and probability predictions on the test set
# # CRITICAL: Include the unique ID (dbgap_subject_id) from the ORIGINAL data (test_m1)
# test_preds <- predict(fit_final_m1, new_data = test_m1, type = "prob") %>%
#   bind_cols(
#     predict(fit_final_m1, new_data = test_m1, type = "class"),
#     # This line ensures the ID and the true outcome are included in the prediction tibble
#     test_m1 %>% dplyr::select(resp_bin, dbgap_subject_id) 
#   )

# Generate class and probability predictions on the test set
test_preds <- predict(fit_final_m1, new_data = test_m1, type = "prob") %>%
  bind_cols(
    predict(fit_final_m1, new_data = test_m1, type = "class"),
    # THIS LINE IS FIXED: Use dbgap_rnaseq_sample
    test_m1 %>% dplyr::select(resp_bin, dbgap_rnaseq_sample) 
  )

# *******************************************************
# ************ CRITICAL STEP: SAVE PREDICTIONS **********
# *******************************************************
test_preds %>%
  # Rename the probability column to clearly identify the model
  rename(prob_Refractory = .pred_Refractory) %>% 
  # Select the essential columns, using the correct ID
  select(dbgap_rnaseq_sample, resp_bin, prob_Refractory) %>%
  # Save to a file specific to the model
  write_csv("predictions_Model1.csv")
# *******************************************************

# The Optimal Threshold Search and subsequent code now follow
# ...

# Optimal Threshold Search
cp <- cutpointr(
  test_preds %>% rename(outcome = resp_bin, probability = .pred_Refractory),
  x = probability, class = outcome, pos_class = "Refractory", neg_class = "CR",
  method = maximize_metric, metric = cutpointr::F1_score 
)
optimal_threshold <- cp$optimal_cutpoint
cat(paste("\nOptimal Threshold (Maximizing F-Measure):", round(optimal_threshold, 4), "\n"))

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

conf_mat(test_preds_optimized, truth = resp_bin, estimate = .pred_class_optimized) %>%
  autoplot(type = "heatmap") 

## 8. ROC & PR curves ------------------------------------------

roc_df <- roc_curve(data  = test_preds, truth = resp_bin, .pred_Refractory)
p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(size = 1) + geom_abline(linetype = "dashed") + coord_equal() +
  labs(title = paste("Model", MODEL_ID, "– ROC curve (Test set)")) + theme_minimal()
print(p_roc) 


## 9. SHAP Feature Importance Analysis (Part 1: Calculation, Bar, and Beeswarm Plot) ------------------------

xgb_engine <- fit_final_m1 %>% extract_fit_engine()

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
# CRITICAL: Pass the feature data as a matrix with NAs.
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

# Create the shapviz object for the standard plot (9c)
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

# --- 9a. SHAP BAR CHART (Mean Absolute SHAP) ---
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


# Extract SHAP data for static ggplot beeswarm
# CRITICAL FIX: Manually convert the raw SHAP matrix to the required long format 
# to avoid the 'sv_data' unexported function error.
shap_data_for_plot <- shap_train %>%
  as.data.frame() %>%
  # Create an index column to pivot around
  mutate(obs_id = row_number()) %>%
  tidyr::pivot_longer(
    cols = -obs_id,
    names_to = "feature",
    values_to = "value"
  ) %>%
  filter(feature %in% top_features) # Filter for the top features

# --- 9b. SHAP Signed Summary Plot (Static Color - NO FEATURE VALUE COLORING, EMPTY CIRCLES) ---
p_shap_signed_static <- ggplot(shap_data_for_plot, aes(x = value, y = feature)) +
  # Add a vertical line at SHAP = 0 (Neutral impact)
  geom_vline(xintercept = 0, color = "gray80", linetype = "dashed") + 
  # Use geom_quasirandom for the beeswarm effect.
  # FIXED: Using shape = 1 for empty circle, color = "black" for the outline.
  ggbeeswarm::geom_quasirandom(shape = 1, color = "black", size = 2.5, alpha = 1, method = "pseudorandom") + 
  scale_y_discrete(limits = rev(top_features)) +
  labs(
    title = paste("Model", MODEL_ID, "– SHAP Signed Summary (Static Color, Empty Circles)"),
    subtitle = "Dot position shows feature impact on prediction (right = Refractory).",
    x = "SHAP Value",
    y = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

cat("\n--- Displaying SHAP Signed Summary (Static Color, Empty Circles) Plot ---\n")
print(p_shap_signed_static) 


# --- 9c. SHAP Signed Summary Plot (Colored by Feature Value - PREVIOUS VERSION) ---
p_shap_signed_colored <- sv_importance(shv, kind = "beeswarm", num_features = N_TOP_FEATURES) +
  labs(
    title = paste("Model", MODEL_ID, "– SHAP Signed Feature Summary (Colored by Value)"),
    subtitle = "Color indicates actual feature value. High values on the right push towards Refractory."
  ) +
  theme_minimal()

cat("\n--- Displaying SHAP Signed Summary (Colored by Value) Plot ---\n")
print(p_shap_signed_colored) 


## 9d. SHAP Distribution by Outcome (CR vs Refractory) with P-values ------------------------

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
  rstatix::wilcox_test(SHAP_Value ~ Outcome) %>% # Explicit call to rstatix
  adjust_pvalue(method = "bonferroni") %>% # Apply multiple testing correction
  rstatix::add_significance("p.adj") %>% # Explicit call to rstatix
  rstatix::add_xy_position(x = "Outcome", dodge = 0.8) %>% # Explicit call to rstatix
  mutate(
    # Use adjusted p-value for significance labeling
    p_text = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ paste("p.adj=", format.pval(p.adj, digits = 2))
    )
  )

cat("\n--- Top", N_TOP_FEATURES, "Features: SHAP P-values (CR vs Refractory) ---\n")
print(p_values_df %>% select(Feature, p, p.adj, p.adj.signif))

# 3. Create the Visualization (Individual Paired Box-and-Jitter Plots)
cat("\n--- Displaying Paired SHAP Distribution Plots ---\n")
# Loop through the top features and generate a separate plot for each
plot_list <- list()
for (feature_name in top_features) {
  
  # Filter data for the current feature
  current_data <- shap_long %>% filter(Feature == feature_name)
  current_p_value <- p_values_df %>% filter(Feature == feature_name)
  
  # Determine max y value for the data points
  max_y_data <- max(current_data$SHAP_Value)
  
  # Calculate 15% of the total range and add it to the max data point for annotation placement
  range_y <- max_y_data - min(current_data$SHAP_Value)
  sig_y_position <- max_y_data + (range_y * 0.15) 
  
  # Update p-value y position for the current plot
  current_p_value <- current_p_value %>% 
    mutate(y.position = sig_y_position)
  
  # Generate the plot
  p <- ggplot(current_data, aes(x = Outcome, y = SHAP_Value, fill = Outcome)) +
    geom_violin(alpha = 0.6, width = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("CR" = "dodgerblue", "Refractory" = "firebrick")) +
    # CRITICAL FIX: Manually expand y-axis limits to ensure annotation is visible
    expand_limits(y = sig_y_position * 1.05) + 
    geom_signif(
      data = current_p_value,
      aes(xmin = xmin, xmax = xmax, annotations = p_text, y = y.position),
      manual = TRUE,
      map_signif_level = FALSE,
      tip_length = 0.01,
      size = 0.5,
      inherit.aes = FALSE # CRITICAL FIX for "object 'Outcome' not found"
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


## 10. Dimensionality Reduction (PCA & UMAP) on Top Features (Targeted Imputation) ---

# 1. Prepare combined data using only the Top Features
X_combined_processed <- bind_rows(X_train_processed, X_test_processed)
X_top_features <- X_combined_processed %>% select(all_of(top_features))

# Get combined outcome data (must match rows of X_top_features)
y_combined_dr <- y_combined

# --- Targeted Imputation for PCA/UMAP ---
cat("\n--- Performing Targeted Imputation for Dimensionality Reduction ---\n")

# Create a function for column-wise median imputation (Package-free solution)
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
# CRITICAL FIX: Corrected typo from data.data.frame to data.frame
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