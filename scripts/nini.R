# ===============================
# inclass-analysis-topN.R
# Compare top 10, 15, 20, 25, 30, 35, 40 proteins
# ===============================

library(tidyverse)
library(randomForest)
library(tidymodels)
library(yardstick)
load('data/biomarker-clean.RData')

# -------------------------------
# Function: top proteins by t-test
# -------------------------------
get_top_ttest <- function(data, n_top = 10){
  ttests_out <- data %>%
    select(-ados) %>%
    pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
    nest(data = c(level, group)) %>%
    mutate(ttest = map(data, ~ t_test(.x, formula = level ~ group, 
                                      order = c('ASD','TD'), 
                                      alternative = 'two-sided'))) %>%
    unnest(ttest) %>%
    arrange(p_value) %>%
    mutate(m = n(),
           hm = log(m) + 1/(2*m) - digamma(1),
           rank = row_number(),
           p.adj = m*hm*p_value/rank)
  
  proteins <- ttests_out %>% slice_min(p.adj, n = n_top) %>% pull(protein)
  return(proteins)
}

# -------------------------------
# Function: top proteins by RF
# -------------------------------
get_top_rf <- function(data, n_top = 10){
  predictors <- data %>% select(-c(group, ados))
  response <- factor(data$group)
  set.seed(101422)
  rf_out <- randomForest(x = predictors, y = response, ntree = 1000, importance = TRUE)
  proteins <- rf_out$importance %>% 
    as_tibble() %>%
    mutate(protein = rownames(rf_out$importance)) %>%
    slice_max(MeanDecreaseGini, n = n_top) %>%
    pull(protein)
  return(proteins)
}

# -------------------------------
# Function: logistic regression & metrics
# -------------------------------
get_logreg_metrics <- function(data, n_top){
  # get top proteins
  proteins_star <- intersect(get_top_ttest(data, n_top), get_top_rf(data, n_top))
  
  # subset data
  biomarker_sstar <- data %>%
    select(group, any_of(proteins_star)) %>%
    mutate(class = (group == 'ASD')) %>%
    select(-group)
  
  # train/test split
  set.seed(101422)
  biomarker_split <- initial_split(biomarker_sstar, prop = 0.8)
  
  # fit logistic regression
  fit <- glm(class ~ ., data = training(biomarker_split), family = 'binomial')
  
  # predict & evaluate
  results <- testing(biomarker_split) %>%
    add_predictions(fit, type = "response") %>%
    mutate(pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE)),
           class = factor(class, levels = c(FALSE, TRUE)))
  
  metrics_res <- metric_set(sensitivity, specificity, accuracy)
  metrics_df <- metrics_res(results, truth = class, estimate = pred_class)
  
  # add N info
  metrics_df$top_proteins <- n_top
  
  return(metrics_df)
}

# -------------------------------
# Run comparison: 10, 15, 20, 25, 30, 35, 40 proteins
# -------------------------------
top_list <- c(10, 15, 20, 25, 30, 35, 40)
all_results <- map_dfr(top_list, ~ get_logreg_metrics(biomarker_clean, .x))

# -------------------------------
# Print comparison table
# -------------------------------
all_results %>%
  select(top_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  arrange(top_proteins)
