library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)

load('../data/biomarker-clean.RData')

# test 2-protein combo with SVM
test_2protein_svm <- function(protein1, protein2, data, seed = 101422) {
  combo_data <- data %>%
    select(group, all_of(c(protein1, protein2))) %>%
    mutate(class = factor(group == 'ASD')) %>%
    select(-group)
  
  set.seed(seed)
  combo_split <- combo_data %>% initial_split(prop = 0.8)
  
  train_combo <- training(combo_split) %>% mutate(class = factor(class))
  fit_svm <- svm(x = train_combo %>% select(-class),
                 y = train_combo$class,
                 kernel = "radial",
                 probability = TRUE)
  
  test_combo <- testing(combo_split) %>% select(-class)
  pred_svm <- predict(fit_svm, test_combo)
  pred_prob <- predict(fit_svm, test_combo, probability = TRUE)
  pred_prob_asd <- attr(pred_prob, "probabilities")[, "TRUE"]
  
  results <- testing(combo_split) %>%
    mutate(pred_class = pred_svm,
           pred_prob = pred_prob_asd,
           truth_class = factor(class)) %>%
    class_metrics(estimate = pred_class,
                  truth = truth_class, 
                  pred_prob,
                  event_level = 'second')
  
  return(results)
}

cat("\n=== RELT + DERM ===\n")
print(test_2protein_svm("RELT", "DERM", biomarker_clean))

cat("\n=== RELT + MRC2 ===\n")
print(test_2protein_svm("RELT", "MRC2", biomarker_clean))
