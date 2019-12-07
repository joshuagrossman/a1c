# analysis.R
# Analyzes cgm, a1c, and demographic data.
# Author: Josh Grossman

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

# Each row is 1 a1c measurement across the 5 datasets. Potentially
# multiple a1c measurements for each patient.
data <- read_csv("data/features.csv")

# 5678 HbA1c, 1296 patients
nrow(data)
length(unique(data$id))

cleaned_df <- load_and_filter(data)

# 3559 HbA1c, 1095 patients
nrow(cleaned_df)
length(unique(cleaned_df$id))

# <10 missing CGM statistics (out of tens of thousands) filled in
# using simple median imputation. Missing ages were aribtrarily cast
# as mean(age) in `load_and_filter`, but will be accounted for later
# on using missingness interaction term per Gelman/Hill.
cleaned_df <- preProcess(cleaned_df, method = "medianImpute") %>%
  predict(cleaned_df) %>% 
  as.data.frame()

train_df <- cleaned_df

######################## MODEL SETUP #################################

covars <- c("mean_bg_full_day", "sd_bg_full_day", 
            "cv_bg_full_day", "percent_very_low_full_day", 
            "percent_low_full_day", 
            "percent_in_target_range_full_day", 
            "percent_in_conservative_target_range_full_day", 
            "percent_high_full_day", "percent_very_high_full_day", 
            "in_target_range_mean_bg_full_day", 
            "high_mean_bg_full_day", 
            "in_target_range_sd_bg_full_day", "high_sd_bg_full_day", 
            "in_target_range_cv_bg_full_day", "high_cv_bg_full_day", 
            "male", "log(age):has_age", "I(log(age)^2):has_age",
            "black", "hispanic")

make_formula <- function(outcome, covars) {
  as.formula(paste0(outcome, " ~ 1 + ", 
                    paste(covars, collapse = " + ")))
}

set_up_model <- function(outcome_name, covars, df) {
  dummied <- dummyVars(make_formula(outcome_name, covars), df,
                       fullRank = TRUE) %>%
    predict(df)

  preProcess(dummied, method = "medianImpute") %>%
    predict(dummied) %>%
    as.data.frame()
}

######### HYPERPARAMETER TUNING, NO PAST HBA1C #######################

pre_processed_regression <- set_up_model("a1c_value", covars, train_df)
outcome_a1c <- train_df$a1c_value

# Ensures that patients are not in both train and validation sets.
tuning_folds <- groupKFold(train_df$id, 5)

tuning_fit_control <- 
  trainControl(method = "repeatedcv",
               index = tuning_folds,
               verboseIter = FALSE,
               allowParallel = TRUE)

# Use LASSO to find predictive coefficients.
lambdas <- 10 ^ seq(-3, 1, length = 100)
tuning_lasso_fit <- train(pre_processed_regression,
                          outcome_a1c,
                          method = "glmnet",
                          tuneGrid = data.frame(alpha = 1,
                                                lambda = lambdas),
                          trControl = tuning_fit_control)
plot(tuning_lasso_fit, xTrans = log, xlab = "log lambda") 
best_lambda_lasso <- tuning_lasso_fit$bestTune$lambda

# # SD of blood glucose and ethnicity look promising.
# coef(lasso_fit$finalModel, lasso_fit$bestTune$lambda)
# 
# # Adding ethnicity does not reduce RMSE much. Note that not a lot of
# # Hispanic patients in the data.
# race_ethnicity_lm_fit <- 
#   train(a1c_value ~ 1 + mean_bg_full_day + black + hispanic,
#         train_df,
#         method = "lm",
#         trControl = fit_control)
# 
# # SD brings down RMSE almost to same level as LASSO.
# race_sd_lm_fit <- 
#   train(a1c_value ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
#         train_df,
#         method = "lm",
#         trControl = fit_control)

# LASSO with all interactions.
tuning_lasso2_fit <- train(a1c_value ~ .*.,
                     bind_cols(a1c_value = outcome_a1c, 
                               pre_processed_regression),
                     method = "glmnet",
                     tuneGrid = data.frame(alpha = 1,
                                           lambda = lambdas),
                     trControl = tuning_fit_control)
plot(tuning_lasso2_fit, xTrans = log, xlab = "log lambda") 
best_lambda_lasso2 <- tuning_lasso2_fit$bestTune$lambda

# Random Forest.
tuning_rf_fit <- train(pre_processed_regression,
                       outcome_a1c,
                       method = "rf",
                       tuneGrid = expand.grid(.mtry=c(1:10)),
                       trControl = tuning_fit_control)
plot(tuning_rf_fit)
best_mtry <- tuning_rf_fit$bestTune$mtry

######### PERFORMANCE METRICS, NO PAST HBA1C  ########################

# Ensures that patients are not in both train and validation sets.
training_folds <- groupKFold(train_df$id, 5)

training_fit_control <- 
  trainControl(method = "repeatedcv",
               index = training_folds,
               verboseIter = FALSE,
               allowParallel = TRUE,
               savePredictions = TRUE)

# GMI is trained on just mean blood glucose, so reproducing.
basic_lm_fit <- train(a1c_value ~ 1 + mean_bg_full_day,
                      train_df,
                      method = "lm",
                      trControl = training_fit_control)

race_lm_fit <- train(a1c_value ~ 1 + mean_bg_full_day + black,
                     train_df,
                     method = "lm",
                     trControl = training_fit_control)

race_sd_lm_fit <- train(a1c_value ~ 1 + mean_bg_full_day + 
                          sd_bg_full_day + black,
                         train_df,
                         method = "lm",
                         trControl = training_fit_control)

lasso_fit <- train(pre_processed_regression,
                   outcome_a1c,
                   method = "glmnet",
                   tuneGrid = data.frame(alpha = 1,
                                         lambda = best_lambda_lasso),
                   trControl = training_fit_control)

lasso2_fit <- train(a1c_value ~ .*.,
                    bind_cols(a1c_value = outcome_a1c, 
                             pre_processed_regression),
                    method = "glmnet",
                    tuneGrid = 
                      data.frame(alpha = 1,
                                 lambda = best_lambda_lasso2),
                    trControl = training_fit_control)

rf_fit <- train(pre_processed_regression,
                       outcome_a1c,
                       method = "rf",
                       tuneGrid = data.frame(mtry = best_mtry),
                       trControl = training_fit_control)

# calculate GMI RMSE on the training folds
rmses <- c()
for (fold in training_folds) {
  gmi_fold <- 3.31 + 0.02392 * train_df$mean_bg_full_day[fold]
  a1c_fold <- outcome_a1c[fold]
  rmse <- sqrt(mean((a1c_fold - gmi_fold)^2))
  rmses <- c(rmses, rmse)
}

gmi_results <- tibble(model = "GMI", 
                      mean_rmse = mean(rmses), 
                      sd_rmse = sd(rmses))

extract_cv_rmse <- function(fit) {
  model <- names(fit)
  rmse <- fit[[model]]$resample$RMSE
  list(model = model, mean_rmse = mean(rmse), sd_rmse = sd(rmse))
}

(performance <- 
    bind_rows(gmi_results, 
              map_dfr(list(
                list(LM_BG = basic_lm_fit), 
                list(LM_BG_RACE = race_lm_fit),
                list(LM_BG_RACE_SD = race_sd_lm_fit),
                list(LASSO = lasso_fit),
                list(LASSO2 = lasso2_fit),
                list(RF = rf_fit)
                ), extract_cv_rmse))) %>% write_csv("lib/analysis/results/full_models_cv.csv")

# Add predictions from best model for later evaluation. Perhaps
# individuals consistently score above or below model, so can use
# residual in future prediction after initial HbA1c lab measurement.
train_df <- mutate(train_df, 
                   a1c_resid = a1c_value - predict(race_sd_lm_fit, 
                                                   train_df))

########## HYPERPARAMETER TUNING, W/ PAST HBA1C ######################

train_df <- extract_last_a1c(train_df)

# Spike of people with at least 70 days between a1c measurements, less
# than that could bias results since a1c is correlated over time.
# Typically a1c measured every 3 months, red blood cells turn over
# every 2-3 months, so this is a reasonable timeframe.
short_train_df <- train_df %>% 
  filter(days_since_last_a1c > 69,
         ! is.na(last_a1c_value))

nrow(short_train_df)
length(unique(short_train_df$id))

short_pre_processed_regression <- 
  set_up_model("a1c_value", c(covars, "last_a1c_value", "last_a1c_resid"), 
               short_train_df)
short_outcome_a1c <- 
  short_train_df$a1c_value

# Ensure that patients are not in both train and validation sets
# during CV.
short_tuning_folds <- groupKFold(short_train_df$id, 5)

short_tuning_fit_control <- 
  trainControl(method = "repeatedcv",
               index = short_tuning_folds,
               verboseIter = FALSE,
               allowParallel = TRUE)

# Use LASSO to find predictive coefficients.
lambdas <- 10 ^ seq(-3, 1, length = 100)
short_tuning_lasso_fit <- train(short_pre_processed_regression,
                          short_outcome_a1c,
                          method = "glmnet",
                          tuneGrid = data.frame(alpha = 1,
                                                lambda = lambdas),
                          trControl = short_tuning_fit_control)
plot(short_tuning_lasso_fit, xTrans = log, xlab = "log lambda") 
short_best_lambda_lasso <- short_tuning_lasso_fit$bestTune$lambda

short_tuning_lasso2_fit <- train(a1c_value ~ .*.,
                           bind_cols(a1c_value = short_outcome_a1c, 
                                     short_pre_processed_regression),
                           method = "glmnet",
                           tuneGrid = data.frame(alpha = 1,
                                                 lambda = lambdas),
                           trControl = short_tuning_fit_control)
plot(short_tuning_lasso2_fit, xTrans = log, xlab = "log lambda") 
short_best_lambda_lasso2 <- short_tuning_lasso2_fit$bestTune$lambda

short_tuning_rf_fit <- train(short_pre_processed_regression,
                       short_outcome_a1c,
                       method = "rf",
                       tuneGrid = expand.grid(.mtry=c(1:10)),
                       trControl = short_tuning_fit_control)
plot(tuning_rf_fit)
short_best_mtry <- short_tuning_rf_fit$bestTune$mtry

########## PERFORMANCE METRICS, W/ PAST HBA1C ########################

short_training_folds <- groupKFold(short_train_df$id, 5)

short_training_fit_control <- 
  trainControl(method = "repeatedcv",
               index = short_training_folds,
               verboseIter = FALSE,
               allowParallel = TRUE,
               savePredictions = TRUE)

short_mrs_fit <-
  train(a1c_value ~ 1 + mean_bg_full_day + black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_training_fit_control)

short_mrs_a1c_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + 
          black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_training_fit_control)

short_mrs_a1c_resid_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + black + 
          sd_bg_full_day + last_a1c_resid,
        short_train_df,
        method = "lm",
        trControl = short_training_fit_control)

short_lasso_fit <- train(short_pre_processed_regression,
                   short_outcome_a1c,
                   method = "glmnet",
                   tuneGrid = data.frame(alpha = 1,
                                         lambda = short_best_lambda_lasso),
                   trControl = short_training_fit_control)

short_lasso2_fit <- train(a1c_value ~ .*.,
                    bind_cols(a1c_value = short_outcome_a1c, 
                              short_pre_processed_regression),
                    method = "glmnet",
                    tuneGrid = 
                      data.frame(alpha = 1,
                                 lambda = short_best_lambda_lasso2),
                    trControl = short_training_fit_control)

short_rf_fit <- train(short_pre_processed_regression,
                short_outcome_a1c,
                method = "rf",
                tuneGrid = data.frame(mtry = short_best_mtry),
                trControl = short_training_fit_control)

short_gmi_rmses <- c()
short_carryover_rmses <- c()

for (fold in short_training_folds) {
  short_gmi <- 3.31 + 0.02392 * short_train_df$mean_bg_full_day[fold]
  short_carryover_a1c <- short_train_df$last_a1c_value[fold]
  short_a1c <- short_train_df$a1c_value[fold]
  
  short_gmi_rmses <- c(short_gmi_rmses, 
                      sqrt(mean((short_gmi - short_a1c)^2)))
  
  short_carryover_rmses <- 
    c(short_carryover_rmses, 
      sqrt(mean((short_carryover_a1c - short_a1c)^2)))
}

short_gmi_results <- tibble(model = "GMI", 
                            mean_rmse = mean(short_gmi_rmses), 
                            sd_rmse = sd(short_gmi_rmses))

short_carryover_results <- tibble(model = "CARRY", 
                            mean_rmse = mean(short_carryover_rmses), 
                            sd_rmse = sd(short_carryover_rmses))

(performance <- 
    bind_rows(short_gmi_results, short_carryover_results) %>% 
    bind_rows(map_dfr(list(
                list(LM_MRS = short_mrs_fit),
                list(`LM_MRS+` = short_mrs_a1c_fit),
                list(`LM_MRS++` = short_mrs_a1c_resid_fit),
                list(LASSO = short_lasso_fit),
                list(LASSO2 = short_lasso2_fit),
                list(RF = short_rf_fit)
              ), extract_cv_rmse))) %>% write_csv("lib/analysis/results/short_models_cv.csv")

############### CROSS-DATASET PERFORMANCE ##############

gmi <- 3.31 + 0.02392 * train_df$mean_bg_full_day
short_gmi <- 3.31 + 0.02392 * short_train_df$mean_bg_full_day

gmi_err <- outcome_a1c - gmi
short_gmi_err <- short_outcome_a1c - short_gmi

sqrt(mean((gmi_err)^2))
sqrt(mean((short_gmi_err)^2))

# Ensure that the same dataset is not in both train and validation
# sets during CV.
cross_fit_control <-
  trainControl(method = "repeatedcv",
               index = groupKFold(train_df$dataset, 5),
               verboseIter = FALSE,
               allowParallel = TRUE)

# rf_cross_fit <- train(pre_processed_regression,
#                 outcome_a1c,
#                 method = "rf",
#                 tuneGrid = data.frame(mtry = best_mtry),
#                 trControl = cross_fit_control)
# 
# cv_rmses <- rf_cross_fit$resample$RMSE
# fold_lengths <- map_int(rf_cross_fit$control$indexOut, length)

race_sd_lm_cross_fit <- train(a1c_value ~ 1 + mean_bg_full_day + 
                          sd_bg_full_day + black,
                        train_df,
                        method = "lm",
                        trControl = cross_fit_control)

cv_rmses <- race_sd_lm_cross_fit$resample$RMSE
fold_lengths <- map_int(race_sd_lm_cross_fit$control$indexOut, length)

# GMI versus MRS+ error across datasets
train_df %>%
  mutate(gmi_err = gmi_err) %>%
  group_by(dataset) %>%
  summarize(n = n(),
            GMI_RMSE = sqrt(mean(gmi_err^2))) %>%
  inner_join(tibble(LM_RMSE = cv_rmses,
                    n = fold_lengths), by = "n") %>%
  mutate(percent_change = (LM_RMSE - GMI_RMSE) / GMI_RMSE) %>% 
  write_csv("lib/analysis/results/full_cdv.csv")


# Ensure that the same dataset is not in both train and validation
# sets during CV.
short_cross_fit_control <- 
  trainControl(method = "repeatedcv",
               index = groupKFold(short_train_df$dataset, 5),
               verboseIter = FALSE,
               allowParallel = TRUE)

mrs_plus_short_cross_fit <- 
  train(a1c_value ~ 1 + last_a1c_value + mean_bg_full_day + 
          black + sd_bg_full_day,
        short_train_df,
        method = "lm",
        trControl = short_cross_fit_control)

cv_rmses <- mrs_plus_short_cross_fit$resample$RMSE
fold_lengths <- map_int(mrs_plus_short_cross_fit$control$indexOut, 
                        length)

# Huge improvements for all the studies, though we have to ignore
# Protocol_F since n is so low, and we can't really trust CMetformin
# for the same reason. This is a strong result that suggests we can do
# better than GMI for sure, especially after a patient has had at
# least one a1c measurement. Perhaps patients now only need to get 1
# a1c per year (instead of the recommended 4).
short_train_df %>% 
  mutate(gmi_err = short_gmi_err) %>% 
  group_by(dataset) %>% 
  summarize(n = n(),
            GMI_RMSE = sqrt(mean(gmi_err^2))) %>% 
  inner_join(tibble(`MRS+_RMSE` = cv_rmses, 
                    n = fold_lengths), by = "n") %>% 
  mutate(percent_change = (`MRS+_RMSE` - GMI_RMSE) / GMI_RMSE) %>% 
  write_csv("lib/analysis/results/short_cdv.csv")

############################## PLOTS ################################

rf_pred <- select(rf_fit$pred, pred, index = rowIndex)

train_df %>% 
  rowid_to_column("index") %>% 
  inner_join(rf_pred, by = "index") %>% 
  mutate(rf_err = a1c_value - pred)