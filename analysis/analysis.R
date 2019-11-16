# analysis.R
# Analyzes cgm, a1c, and demographic data.
# Author: Josh Grossman

################################################################################
#
# IMPORTANT NOTE TO READER: 
# The full reposistory of code used for this project can be found at the 
# following link: https://github.com/joshuagrossman/a1c
#
# The code attached to the project is only the code required for analysis.
# There are *many* other files written to clean and process the combined
# datasets used in this project (several GBs).
#
# The code to create the best regression and classification model has been 
# surrounded with 3 lines of hashes above and below. 
#
################################################################################

source("lib/analysis/helpers.R")

library(caret)
library(glmnet)
library(lubridate)
library(randomForest)
library(furrr)
library(e1071)
library(ROCR)

# For fitting random forest models quickly.
library(doParallel)
cl <- makePSOCKcluster(availableCores() %/% 2 + 1)
registerDoParallel(cl)

set.seed(1)

TRAIN_FRAC <- 0.8

# Each row is 1 a1c measurement across the 5 datasets. Potentially multiple a1c
# measurements for each patient.
data <- read_csv("data/features.csv")

cleaned_df <- load_and_filter(data)

# <10 missing CGM statistics (out of tens of thousands) filled in using simple
# median imputation.
# Missing ages were cast as mean(age) in `load_and_filter`, and will be 
# accounted for later on using interaction term per Gelman/Hill chapter.
cleaned_df <- preProcess(cleaned_df, method = "medianImpute") %>%
  predict(cleaned_df) %>% 
  as.data.frame()

# Generate train and test data.
# Don't want the same patient to exist in the train and test data because
# multiple outcomes per patient.
training_ids <- 
  cleaned_df %>% 
  select(id) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  sample_frac(TRAIN_FRAC)

training <-
  cleaned_df %>%
  rowid_to_column("index") %>% 
  right_join(training_ids, by = "id") %>% 
  pull("index")

train_df <- cleaned_df[training, ]
test_df <- cleaned_df[-training, ]

# # Check that there are no matching IDs in train and test sets.
# nrow(cleaned_df)
# nrow(train_df) + nrow(test_df)
# nrow(inner_join(train_df, test_df, by = "id"))

# write_csv(train_df, "data/train.csv")
# write_csv(test_df, "data/test.csv")

#train_df <- read_csv("data/train.csv")
#test_df <- read_csv("data/test.csv")

# Baseline models for regression and classification. Current "gold-standard".
gmi <- 3.31 + 0.02392 * train_df$mean_bg_full_day
gmi_unhealthy <- as.numeric(gmi > 7.5)

# Using RMSE for regression evaluation.
outcome_a1c <- train_df$a1c_value
gmi_rmse <- sqrt(mean((outcome_a1c - gmi)^2))

# Using AUC for classification evaluation.
unhealthy_a1c <- train_df$unhealthy_a1c
gmi_auc <- performance(prediction(gmi_unhealthy, unhealthy_a1c), "auc")@y.values[[1]]

######################## MODEL SETUP ###########################################

train_df$unhealthy_a1c <- if_else(train_df$unhealthy_a1c == 1, 
                                  "Unhealthy", "Healthy")

covars <- c("mean_bg_full_day", "sd_bg_full_day", 
            "cv_bg_full_day", "percent_very_low_full_day", 
            "percent_low_full_day", "percent_in_target_range_full_day", 
            "percent_in_conservative_target_range_full_day", 
            "percent_high_full_day", "percent_very_high_full_day", 
            "in_target_range_mean_bg_full_day", "high_mean_bg_full_day", 
            "in_target_range_sd_bg_full_day", "high_sd_bg_full_day", 
            "in_target_range_cv_bg_full_day", "high_cv_bg_full_day", 
            "male", "log(age):has_age", "I(log(age)^2):has_age",
            "black", "hispanic")

make_formula <- function(outcome, covars) {
  as.formula(paste0(outcome, " ~ 1 + ", paste(covars, collapse = " + ")))
}

set_up_model <- function(outcome_name, covars) {
  dummied <- dummyVars(make_formula(outcome_name, covars), train_df,
                       fullRank = TRUE) %>%
    predict(train_df)

  preProcess(dummied, method = "medianImpute") %>%
    predict(dummied) %>%
    as.data.frame()
}

pre_processed_regression <- set_up_model("a1c_value", covars)
pre_processed_classification <- set_up_model("unhealthy_a1c", covars)

######################## REGRESSION: PART 1 ####################################

# Ensures that patients are not in both train and validation sets.
fit_control <- trainControl(method = "repeatedcv",
                            index = groupKFold(train_df$id, 5),
                            verboseIter = FALSE,
                            allowParallel = TRUE)

# Sequential model fitting.

# GMI is trained on just mean blood glucose, so reproducing.
basic_lm_fit <- train(a1c_value ~ 1 + mean_bg_full_day,
                      train_df,
                      method = "lm",
                      trControl = fit_control)

basic_lm_fit$finalModel

# There is evidence that a1c differs across race, so adding into model.
race_lm_fit <- train(a1c_value ~ 1 + mean_bg_full_day + black,
                     train_df,
                     method = "lm",
                     trControl = fit_control)

# Using LASSO to find predictive coefficients.
lambdas <- 10 ^ seq(-3, 1, length = 100)
lasso_fit <- train(pre_processed_regression,
                   outcome_a1c,
                   method = "glmnet",
                   tuneGrid = data.frame(alpha = 1,
                                         lambda = lambdas),
                   trControl = fit_control)
plot(lasso_fit, xTrans = log) 

# SD of blood glucose and ethnicity look promising. Other covariates are 
# almost too "specific", so will try to reduce RMSE to LASSO level by adding 
# the most interpretable covariates as possible.
coef(lasso_fit$finalModel, lasso_fit$bestTune$lambda)

# Adding ethnicity does not reduce RMSE much. Note that not a lot of Hispanic
# patients in the data.
race_ethnicity_lm_fit <- 
  train(a1c_value ~ 1 + mean_bg_full_day + black + hispanic,
        train_df,
        method = "lm",
        trControl = fit_control)

# SD brings down RMSE almost to same level as LASSO.
race_sd_lm_fit <- 
  train(a1c_value ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
        train_df,
        method = "lm",
        trControl = fit_control)

# LASSO with all interactions does not improve much on LASSO.
lasso2_fit <- train(a1c_value ~ .*.,
                    bind_cols(a1c_value = outcome_a1c, pre_processed_regression),
                    method = "glmnet",
                    tuneGrid = data.frame(alpha = 1,
                                          lambda = lambdas),
                    trControl = fit_control)
plot(lasso2_fit, xTrans = log) 

# Random Forest does not improve much beyond LASSO.
rf_fit <- train(pre_processed_regression,
                outcome_a1c,
                method = "rf",
                tuneGrid = expand.grid(.mtry=c(1:15)),
                trControl = fit_control)
plot(rf_fit)

extract_cv_rmse <- function(fit) {
  mean(fit$resample$RMSE)
}

# LASSO performs just as well as LASSO^2 and RF. 
# LM w/ race is pretty good though, and more "interpretable".
# LM w/ race and SD is also interpretable, and gets close to LASSO. Best model.
list(gmi_rmse = gmi_rmse) %>% 
  append(
    map(list(basic_lm_cv_rmse = basic_lm_fit, 
             race_lm_cv_rmse = race_lm_fit,
             race_ethnicity_lm_cv_rmse = race_ethnicity_lm_fit,
             race_sd_lm_cv_rmse = race_sd_lm_fit,
             lasso_cv_rmse = lasso_fit,
             lasso2_cv_rmse = lasso2_fit,
             rf_cv_rmse = rf_fit),
             extract_cv_rmse))

# Add predictions from best model for later evaluation. Perhaps individuals
# consistently score above or below model, so can use residual in future 
# prediction after initial HbA1c lab measurement.
train_df <- mutate(train_df, 
                   a1c_resid = a1c_value - predict(race_sd_lm_fit, train_df))

######################## CLASSIFICATION ########################################

# Ensures that patients are not in both train and validation sets.
class_fit_control <- trainControl(method = "repeatedcv",
                            index = groupKFold(train_df$id, 5),
                            verboseIter = FALSE,
                            allowParallel = TRUE,
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary)

basic_glm_fit <- train(unhealthy_a1c ~ 1 + mean_bg_full_day,
                      train_df,
                      method = "glm",
                      family = binomial(),
                      metric = "ROC",
                      trControl = class_fit_control)

race_glm_fit <- train(unhealthy_a1c ~ 1 + mean_bg_full_day + black,
                     train_df,
                     method = "glm",
                     family = binomial(),
                     metric = "ROC",
                     trControl = class_fit_control)

lambdas <- 10 ^ seq(-3, 1, length = 100)

class_lasso_fit <- train(pre_processed_classification,
                   train_df$unhealthy_a1c,
                   method = "glmnet",
                   family = "binomial",
                   tuneGrid = data.frame(alpha = 1,
                                         lambda = lambdas),
                   metric = "ROC",
                   trControl = class_fit_control)
plot(class_lasso_fit, xTrans = log) 

coef(class_lasso_fit$finalModel, class_lasso_fit$bestTune$lambda)

race_ethnicity_glm_fit <- 
  train(unhealthy_a1c ~ 1 + mean_bg_full_day + black + hispanic,
        train_df,
        method = "glm",
        family = binomial(),
        metric = "ROC",
        trControl = class_fit_control)

################################################################################
################################################################################
################################################################################

# BEST CLASSIFICATION MODEL (logistic MRS model)

race_sd_glm_fit <- 
  train(unhealthy_a1c ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
        train_df,
        method = "glm",
        metric = "ROC",
        trControl = class_fit_control)

################################################################################
################################################################################
################################################################################

class_lasso2_fit <- train(unhealthy_a1c ~ .*.,
                    bind_cols(unhealthy_a1c = train_df$unhealthy_a1c, pre_processed_classification),
                    method = "glmnet",
                    family = "binomial",
                    tuneGrid = data.frame(alpha = 1,
                                          lambda = lambdas),
                    metric = "ROC",
                    trControl = class_fit_control)
plot(class_lasso2_fit, xTrans = log) 

class_rf_fit <- train(pre_processed_classification,
                train_df$unhealthy_a1c,
                method = "rf",
                family = "binomial",
                tuneGrid = expand.grid(.mtry=c(1:15)),
                metric = "ROC",
                trControl = class_fit_control)
plot(class_rf_fit)

extract_cv_auc <- function(fit) {
  mean(fit$resample$ROC)
}

# Nothing is *that* much better than GMI in terms of accuracy.
# Again, though, a linear model with race/mean BG/sd BG works just as well
# as the more complicated models, so will call this the best model.
list(gmi_auc = gmi_auc) %>% 
  append(
    map(list(basic_glm_cv_auc = basic_glm_fit, 
             race_glm_cv_auc = race_glm_fit,
             race_ethnicity_glm_cv_auc = race_ethnicity_glm_fit,
             race_sd_glm_cv_auc = race_sd_glm_fit,
             class_lasso_cv_auc = class_lasso_fit,
             class_lasso2_cv_auc = class_lasso2_fit,
             class_rf_cv_auc = class_rf_fit),
        extract_cv_auc))

# Adding residuals for potential later modeling.
train_df <- mutate(train_df, 
                   unhealthy_resid = as.integer(unhealthy_a1c == predict(race_sd_glm_fit, train_df)))

######################## REGRESSION: PART 2 ####################################

train_df <- extract_last_a1c(train_df)

# Spike of people with at least 70 days between a1c measurements, less than
# that could skew results since a1c is correlated over time.
# Typically a1c measured every 3 months, red blood cells turn over every 2-3
# months, so this is a reasonable timeframe.
short_train_df <- train_df %>% 
  filter(days_since_last_a1c > 69,
         ! is.na(last_a1c_value))

# Ensure that patients are not in both train and validation sets during CV.
short_fit_control <- trainControl(method = "repeatedcv",
                                  index = groupKFold(short_train_df$id, 5),
                                  verboseIter = FALSE,
                                  allowParallel = TRUE)

# Baseline models are GMI and simply predicting the same a1c as last time.
short_gmi <- 3.31 + 0.02392 * short_train_df$mean_bg_full_day
short_carryover_a1c <- short_train_df$last_a1c_value
short_outcome_a1c <- short_train_df$a1c_value

short_gmi_rmse <- sqrt(mean((short_outcome_a1c - short_gmi)^2))
carryover_rmse <- sqrt(mean((short_outcome_a1c - short_carryover_a1c)^2))

# MRS model.
former_fit <-
  train(a1c_value ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_fit_control)

################################################################################
################################################################################
################################################################################

# BEST REGRESSION MODEL (MRS+)

last_a1c_former_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_fit_control)

################################################################################
################################################################################
################################################################################

# Best model augmented with the previous a1c value and the residual from
# previous prediction. (MRS++)
resid_last_a1c_former_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + black + sd_bg_full_day + last_a1c_resid,
        short_train_df,
        method = "lm",
        trControl = short_fit_control)

# Including the previous a1c substantially improves the RMSE. Adding the
# residual doesn't do much beyond that.
list(short_gmi_rmse = short_gmi_rmse,
     carryover_rmse = carryover_rmse) %>% 
  append(
    map(list(former_cv_rmse = former_fit, 
             last_a1c_former_cv_rmse = last_a1c_former_fit,
             resid_last_a1c_former_cv_rmse = resid_last_a1c_former_fit),
        extract_cv_rmse))

#################### CROSS-DATASET VALIDATION: REGRESSION ######################

# Ensure that the same dataset is not in both train and validation sets during CV.
short_cross_fit_control <- trainControl(method = "repeatedcv",
                                  index = groupKFold(short_train_df$dataset, 5),
                                  verboseIter = FALSE,
                                  allowParallel = TRUE)

mrs_plus_short_cross_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_cross_fit_control)

cv_rmses <- mrs_plus_short_cross_fit$resample$RMSE
fold_lengths <- map_int(mrs_plus_short_cross_fit$control$indexOut, length)

# Huge improvements for all the studies, though we have to ignore Protocol_F
# since n is so low, and we can't really trust CMetformin for the same reason.
# This is a strong result that suggests we can do better than GMI for sure,
# especially after a patient has had at least one a1c measurement. Perhaps
# patients now only need to get 1 a1c per year (instead of the recommended 4).
short_train_df %>% 
  mutate(gmi_err = a1c_value - (3.31 + 0.02392 * mean_bg_full_day)) %>% 
  group_by(dataset) %>% 
  summarize(n = n(),
            GMI_RMSE = sqrt(mean(gmi_err^2))) %>% 
  inner_join(tibble(`MRS+_RMSE` = cv_rmses, n = fold_lengths), by = "n") %>% 
  mutate(percent_change = (`MRS+_RMSE` - GMI_RMSE) / GMI_RMSE)

#################### CROSS-DATASET VALIDATION: CLASSIFICATION ##################

# Ensure that the same dataset is not in both train and validation sets during CV.
cross_fit_control <- trainControl(method = "repeatedcv",
                                  index = groupKFold(train_df$dataset, 5),
                                  verboseIter = FALSE,
                                  allowParallel = TRUE,
                                  classProbs = TRUE,
                                  summaryFunction = twoClassSummary)

mrs_class_cross_fit <- 
  train(unhealthy_a1c ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
        train_df,
        method = "glm",
        family = binomial(),
        metric = "ROC",
        trControl = cross_fit_control)

cv_auc <- mrs_class_cross_fit$resample$ROC
cv_auc

#################### ROC CURVE FOR LOGISTIC MRS MODEL ##########################

set.seed(1)
roc_train <- sample(1:nrow(train_df), round(0.8 * nrow(train_df)))
roc_train_df <- train_df[roc_train,]
roc_validation_df <- train_df[-roc_train,]

roc_fit <- glm(factor(unhealthy_a1c) ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
               roc_train_df, family = binomial())
pred <- predict(roc_fit, newdata = roc_validation_df, type = "response")
true <- as.numeric(factor(roc_validation_df$unhealthy_a1c)) - 1

roc_curve <- performance(prediction(pred, true),"fpr","fnr")

roc_df <- tibble(fnr = unlist(roc_curve@x.values), 
                 fpr = unlist(roc_curve@y.values)) %>% 
  mutate(fnr_fpr_ratio = fnr/fpr,
         dist_from_10 = abs(fnr_fpr_ratio - 10)) %>% 
  arrange(dist_from_10)

head(roc_df)

###################### DEPRECATED ##############################################

# ### Back to classification
# 
# # Baseline models are GMI and using the same classification as the last a1c.
# short_gmi_unhealthy <- factor(as.integer(short_gmi < 7.5))
# short_carryover_unhealthy <- factor(as.integer(short_carryover_a1c < 7.5))
# short_unhealthy_a1c <- short_train_df$unhealthy_a1c
# 
# short_gmi_acc <- mean(short_unhealthy_a1c == short_gmi_unhealthy)
# carryover_acc <- mean(short_unhealthy_a1c == short_carryover_unhealthy)
# 
# # Best model from before.
# class_former_fit <-
#   train(unhealthy_a1c ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
#         short_train_df,
#         method = "glm",
#         family = binomial(),
#         trControl = short_fit_control)
# 
# # Best model augmented with the previous a1c value.
# last_a1c_class_former_fit <- 
#   train(unhealthy_a1c ~ 1 + last_unhealthy_a1c + mean_bg_full_day + 
#           black + sd_bg_full_day,
#         short_train_df,
#         method = "glm",
#         family = binomial(),
#         trControl = short_fit_control)
# 
# # Best model augmented with the previous a1c value and the residual from
# # previous prediction.
# resid_last_a1c_class_former_fit <- 
#   train(unhealthy_a1c ~ 1 + last_unhealthy_a1c + mean_bg_full_day + black + 
#           sd_bg_full_day + last_unhealthy_resid,
#         short_train_df,
#         method = "glm",
#         family = binomial(),
#         trControl = short_fit_control)
# 
# # Including the prior a1c is a solid improvement over GMI.
# # Interesting, GMI is no worse than the best model that does not include the
# # prior a1c.
# list(short_gmi_acc = short_gmi_acc,
#      carryover_acc = carryover_acc) %>% 
#   append(
#     map(list(former_cv_acc = class_former_fit, 
#              last_a1c_former_cv_acc = last_a1c_class_former_fit,
#              resid_last_a1c_former_cv_acc = resid_last_a1c_class_former_fit),
#         extract_cv_accuracy))


# fold_lengths <- map_int(mrs_class_cross_fit$control$indexOut, length)
# 
# unhealthy_a1c_raw <- as.numeric(train_df$unhealthy_a1c) - 1
# 
# train_df %>% 
#   mutate(gmi = 3.31 + 0.02392 * mean_bg_full_day,
#          GMI_AUC = abs(unhealthy_a1c_raw - (gmi < 7.5))) %>% 
#   group_by(dataset) %>% 
#   summarize(n = n(),
#             GMI_ACC = 1 - mean(gmi_err01)) %>% 
#   inner_join(tibble(`MRS_ACC` = cv_acc, n = fold_lengths), by = "n") %>% 
#   mutate(percent_change = (`MRS_ACC` - GMI_ACC) / GMI_ACC)