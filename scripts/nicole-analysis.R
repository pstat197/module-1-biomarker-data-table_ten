library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)
load('data/biomarker-clean.RData')

predictors <- biomarker_clean %>%
  select(-c(group, ados)) %>%
  as.matrix()

response <- if_else(biomarker_clean$group == "ASD", 1, 0)

# partition into training and test set
set.seed(101422)
split <- initial_split(biomarker_clean, prop = 0.8)
train_data <- training(split)
test_data  <- testing(split)

x_train <- train_data %>% select(-c(group, ados)) %>% as.matrix()
y_train <- if_else(train_data$group == "ASD", 1, 0)

x_test <- test_data %>% select(-c(group, ados)) %>% as.matrix()
y_test <- if_else(test_data$group == "ASD", 1, 0)

#fit lasso regularization
set.seed(101422)
cvfit <- cv.glmnet(
  x_train, y_train,
  family = "binomial",
  alpha = 1,       # LASSO
  nfolds = 5
)

# lambda that minimizes cross-validated error
lambda_min <- cvfit$lambda.min

#importance
lasso_coefs <- coef(cvfit, s = "lambda.min")
proteins_lasso <- rownames(lasso_coefs)[lasso_coefs[,1] != 0]
proteins_lasso <- setdiff(proteins_lasso, "(Intercept)")

proteins_lasso

# predicted probabilities
test_data_lasso <- test_data %>%
  mutate(
    pred = as.numeric(predict(cvfit, newx = x_test, s = lambda_min, type = "response")),
    preds = factor(if_else(pred > 0.3, "ASD", "TD"), levels = c("TD", "ASD")),
    classes = factor(if_else(y_test == 1, "ASD", "TD"), levels = c("TD", "ASD"))
  )

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

class_metrics(
  test_data_lasso,
  estimate = preds,
  truth = classes,
  pred,
  event_level = "second"
)

