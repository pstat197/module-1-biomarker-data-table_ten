# install.packages("randomForest")
library(tidyverse)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load('../data/biomarker-clean.RData')
set.seed(101422)

# ============================================================================
# BASELINE (from in-class)
# ============================================================================

# T-tests
ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, ~ t_test(.x, formula = level ~ group,
                                    order = c('ASD', 'TD'),
                                    alternative = 'two-sided',
                                    var.equal = F))) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(), p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>% slice_min(p.adj, n = 10) %>% pull(protein)

# Random Forest
predictors <- biomarker_clean %>% select(-c(group, ados))
response <- biomarker_clean %>% pull(group) %>% factor()
rf_out <- randomForest(x = predictors, y = response, ntree = 1000, importance = T)
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

# Baseline: intersection
proteins_baseline <- intersect(proteins_s1, proteins_s2)

# Evaluation function
evaluate_panel <- function(protein_set) {
  data <- biomarker_clean %>%
    select(group, any_of(protein_set)) %>%
    mutate(class = (group == 'ASD')) %>%
    select(-group)
  
  split <- initial_split(data, prop = 0.8)
  fit <- glm(class ~ ., data = training(split), family = 'binomial')
  
  metrics <- testing(split) %>%
    add_predictions(fit, type = 'response') %>%
    metric_set(accuracy, sensitivity, specificity, roc_auc)(
      estimate = factor(pred > 0.5),
      truth = factor(class), pred,
      event_level = 'second')
  
  return(metrics)
}

baseline_metrics <- evaluate_panel(proteins_baseline)
cat("\n=== BASELINE ===\n")
cat("Proteins:", length(proteins_baseline), "\n")
print(baseline_metrics)

# ============================================================================
# SIMPLIFIED PANEL: Top K by Combined Ranking
# ============================================================================

# Combine importance from both methods
rf_importance <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance),
         rf_rank = rank(-MeanDecreaseGini))

combined_rank <- ttests_out %>%
  select(protein, p.adj) %>%
  mutate(ttest_rank = rank(p.adj)) %>%
  left_join(rf_importance, by = "protein") %>%
  mutate(avg_rank = (rf_rank + ttest_rank) / 2) %>%
  arrange(avg_rank)

# Test different panel sizes
cat("\n=== TESTING SIMPLIFIED PANELS ===\n")
results <- tibble()

for(k in 2:8) {
  proteins <- combined_rank %>% slice_head(n = k) %>% pull(protein)
  metrics <- evaluate_panel(proteins)
  
  acc <- metrics %>% filter(.metric == "accuracy") %>% pull(.estimate)
  auc <- metrics %>% filter(.metric == "roc_auc") %>% pull(.estimate)
  
  results <- bind_rows(results, tibble(n_proteins = k, accuracy = acc, roc_auc = auc))
  cat("K =", k, "| Accuracy:", round(acc, 3), "| AUC:", round(auc, 3), "\n")
}

# Find best simplified panel (smallest that's within 5% of baseline)
baseline_acc <- baseline_metrics %>% filter(.metric == "accuracy") %>% pull(.estimate)
baseline_auc <- baseline_metrics %>% filter(.metric == "roc_auc") %>% pull(.estimate)

best <- results %>%
  filter(accuracy >= baseline_acc - 0.05, 
         roc_auc >= baseline_auc - 0.05) %>%
  slice_min(n_proteins, n = 1)

cat("\n=== RECOMMENDED SIMPLIFIED PANEL ===\n")
cat("Size:", best$n_proteins, "proteins (vs", length(proteins_baseline), "baseline)\n")
cat("Accuracy:", round(best$accuracy, 3), "(vs", round(baseline_acc, 3), "baseline)\n")
cat("AUC:", round(best$roc_auc, 3), "(vs", round(baseline_auc, 3), "baseline)\n")

recommended_proteins <- combined_rank %>% slice_head(n = best$n_proteins) %>% pull(protein)
cat("\nProteins:", paste(recommended_proteins, collapse = ", "), "\n")

# Quick visualization
ggplot(results, aes(x = n_proteins, y = accuracy)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = baseline_acc, linetype = "dashed", color = "red") +
  labs(title = "Simplified Panel Performance",
       x = "Number of Proteins", y = "Accuracy") +
  theme_minimal()

ggsave("results/simplified-panel-plot.png", width = 6, height = 4)