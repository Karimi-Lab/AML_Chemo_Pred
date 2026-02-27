## ============================================================
## MODEL COMPARISON SCRIPT (Model 1, 2, 3)
## Generates Overlaid ROC and PR Curves and a summary table.
##
## PREREQUISITE: Must have the following files in your working directory:
## 1. predictions_Model1.csv
## 2. predictions_Model2.csv
## 3. predictions_Model3.csv
## ============================================================

# 0. Libraries ------------------------------------------------
# Ensure these are installed: install.packages(c("tidymodels", "dplyr", "readr", "ggplot2", "patchwork", "tibble"))
library(tidymodels)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork) # For arranging plots
library(tibble)    # For easy data frame manipulation
tidymodels_prefer()

# 1. Load and combine prediction data -------------------------

# Define column types explicitly, resolving the namespace conflict for col_factor
col_spec <- readr::cols(
  dbgap_rnaseq_sample = readr::col_character(),
  resp_bin = readr::col_factor(levels = c("CR", "Refractory")), # RESOLVED CONFLICT
  prob_Refractory = readr::col_double()
)

cat("Loading predictions...\n")

# Load predictions for each model. 
# We rename any 'dbgap_...' column to the standard 'dbgap_rnaseq_sample' for binding consistency.
tryCatch({
  preds_m1 <- readr::read_csv("predictions_Model1.csv", col_types = col_spec) %>% 
    rename_with(~"dbgap_rnaseq_sample", dplyr::contains("dbgap_")) %>%
    mutate(Model = "Model 1: Baseline")
  
  preds_m2 <- readr::read_csv("predictions_Model2.csv", col_types = col_spec) %>% 
    rename_with(~"dbgap_rnaseq_sample", dplyr::contains("dbgap_")) %>%
    mutate(Model = "Model 2: RNA-seq Only")
  
  preds_m3 <- readr::read_csv("predictions_Model3.csv", col_types = col_spec) %>% 
    rename_with(~"dbgap_rnaseq_sample", dplyr::contains("dbgap_")) %>%
    mutate(Model = "Model 3: RNA-seq + TEF")
  
}, error = function(e) {
  cat("Error loading files. Please ensure all three files exist and contain columns like 'resp_bin', 'prob_Refractory', and one unique ID column.\n")
  stop(e)
})


# Combine all predictions into one dataframe
all_predictions <- dplyr::bind_rows(preds_m1, preds_m2, preds_m3) %>%
  # Ensure the outcome is correctly leveled
  mutate(resp_bin = factor(resp_bin, levels = c("CR", "Refractory")))

cat("Combined data loaded. Total predictions:", nrow(all_predictions), "\n")


# 2. Calculate Scalar Metrics (ROC AUC and PR AUC) and create Plot Labels ------------

cat("\n=== Comparative Scalar Metrics (Test Set) ===\n")
metric_set_comp <- metric_set(roc_auc, pr_auc)

all_metrics_df <- all_predictions %>%
  group_by(Model) %>%
  metric_set_comp(
    truth = resp_bin,
    prob_Refractory, 
    event_level = "second"
  ) %>%
  ungroup() %>%
  # Pivot wider to easily combine the AUC values into a label
  tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
  # Create a descriptive label for the plot legend
  mutate(
    Model_Label = paste0(
      Model,
      " (ROC AUC: ", format(round(roc_auc, 3), nsmall = 3),
      "; PR AUC: ", format(round(pr_auc, 3), nsmall = 3),
      ")"
    )
  )

# Print the final summary table
metrics_table <- all_metrics_df %>% 
  select(Model, `ROC AUC` = roc_auc, `PR AUC` = pr_auc) %>% 
  arrange(desc(`ROC AUC`)) %>%
  # Use tibble to display cleanly
  as_tibble()

cat("\n")
print(metrics_table)
cat("\n")

# Store the Model_Label for consistent plotting
model_labels <- all_metrics_df %>% 
  select(Model, Model_Label) %>% 
  distinct()


# 3. Generate Overlaid ROC Plot -------------------------------

# Calculate ROC curves for all models
roc_curves_df <- all_predictions %>%
  group_by(Model) %>%
  roc_curve(truth = resp_bin, prob_Refractory) %>%
  ungroup() %>%
  left_join(model_labels, by = "Model") # Join with custom labels

# Define colors for the plot (using a color-blind friendly palette)
model_colors <- RColorBrewer::brewer.pal(3, "Dark2")

p_roc_overlap <- ggplot(roc_curves_df, aes(x = 1 - specificity, y = sensitivity, color = Model_Label)) +
  geom_path(size = 1.2) +
  geom_abline(lty = 2, alpha = 0.5, color = "gray50") + # Diagonal reference line
  coord_equal() +
  scale_color_manual(values = model_colors) +
  labs(
    title = "A. Comparative Receiver Operating Characteristic (ROC) Curves",
    subtitle = "Performance of all models on the held-out Test Set",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Model Performance"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

cat("--- Displaying Overlaid ROC Curve ---\n")
print(p_roc_overlap) 


# 4. Generate Overlaid PR Plot --------------------------------

# Calculate PR curves for all models
pr_curves_df <- all_predictions %>%
  group_by(Model) %>%
  pr_curve(truth = resp_bin, prob_Refractory) %>%
  ungroup() %>%
  left_join(model_labels, by = "Model") # Join with custom labels

p_pr_overlap <- ggplot(pr_curves_df, aes(x = recall, y = precision, color = Model_Label)) +
  geom_path(size = 1.2) +
  scale_color_manual(values = model_colors) +
  labs(
    title = "B. Comparative Precision-Recall (PR) Curves",
    x = "Recall",
    y = "Precision",
    color = "Model Performance"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

cat("\n--- Displaying Overlaid PR Curve ---\n")
print(p_pr_overlap) 


# 5. Final Combined Plot (Optional) ---------------------------

# Use patchwork to place plots side-by-side with a shared legend
p_combined <- p_roc_overlap + p_pr_overlap + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Comparison of Predictive Model Performance (Test Set)",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  ) & theme(legend.position = 'bottom')

cat("\n--- Displaying Combined ROC and PR Plot ---\n")
print(p_combined)