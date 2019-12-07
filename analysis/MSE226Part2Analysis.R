# MSE226Part2Analsis.R
# Analyzes cgm, a1c, and demographic data for part 2.
# Author: Josh Grossman

######################################################################
#
# IMPORTANT NOTE TO READER: The full reposistory of code used for this
# project can be found at the following link:
# https://github.com/joshuagrossman/a1c
#
# The code attached to the project is only the code required for
# analysis. There are *many* other files written to clean and process
# the combined datasets used in this project (several GBs of data).
#
# Most of the datasets were obtained from 
# https://public.jaeb.org/t1dx/stdy, though the FLEX data was
# obtained through collaboration with a group at UNC.
#
######################################################################

source("lib/analysis/MSE226Part1Helpers.R")

library(tidyverse)

train_df <- 
  read_csv("data/train.csv") %>% 
  extract_last_a1c() %>% 
  filter(! is.na(last_a1c_value))

test_df <- read_csv("data/test.csv") %>% 
  extract_last_a1c() %>% 
  filter(! is.na(last_a1c_value))

######################################################################

#####
# 1 #
#####

# Calculating test error for the MRS+ model.
# Estimated test error from CV was 0.47 in part 1.

mrs_plus_formula <- a1c_value ~ 1 + last_a1c_value + 
                    mean_bg_full_day + black + sd_bg_full_day

train_fit <- lm(mrs_plus_formula, train_df)

preds <- predict(train_fit, newdata = test_df)

sqrt(mean((test_df$a1c_value - preds)^2))

# Test error is 0.45. 
# The estimated test error was a good approximation.

# Calculating FPR and FNR where FPR * 10 = FNR.
# These values were estimated at 35.9% FNR and 3.6% FPR.

roc_fit <- glm(factor(unhealthy_a1c) ~ 1 + mean_bg_full_day + 
                 black + sd_bg_full_day,
               train_df, family = binomial())
pred <- predict(roc_fit, newdata = test_df, 
                type = "response")
true <- as.numeric(factor(test_df$unhealthy_a1c)) - 1

roc_curve <- performance(prediction(pred, true),"fpr","fnr")

roc_df <- tibble(fnr = unlist(roc_curve@x.values), 
                 fpr = unlist(roc_curve@y.values)) %>% 
  mutate(fnr_fpr_ratio = fnr/fpr,
         dist_from_10 = abs(fnr_fpr_ratio - 10)) %>% 
  arrange(dist_from_10)

head(roc_df)

# On test set, 51.9% FNR and 5.1% FPR.
# The estimates of these values were overconfident.

######################################################################
  
##### 
# 2 #
#####

#####
# a #
#####

summary(train_fit)

#####
# b #
#####

test_fit <- lm(mrs_plus_formula, test_df)

summary(test_fit)

inner_join(broom::tidy(train_fit), broom::tidy(test_fit),
           by = "term") %>% 
  map_if(is.numeric, ~ round(.x, 3)) %>% 
  as_tibble() %>% 
  mutate(
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "last_a1c_value" ~ "Prior HbA1c",
      term == "mean_bg_full_day" ~ "Mean BG",
      term == "black" ~ "Black",
      term == "sd_bg_full_day" ~ "SD BG",
      TRUE ~ NA_character_
    )
  ) %>% 
  rename(
    Coefficient = term,
    `Tr est` = estimate.x,
    `Tr SE` = std.error.x,
    `Tr t` = statistic.x,
    `Tr p` = p.value.x,
    `Tst est` = estimate.y,
    `Tst SE` = std.error.y,
    `Tst t` = statistic.y,
    `Tst p` = p.value.y
  ) %>% 
  write_csv("lib/analysis/inference-results-ab.csv")

#####
# c #
#####

set.seed(1)

calc_boot_coefs <- function(df) {
  boot_df <- sample_frac(df, replace = TRUE)
  boot_fit <- lm(mrs_plus_formula, boot_df)
  coef(boot_fit)
}

nboot <- 5000

boot_coefs <- 
  map(1:nboot, ~ calc_boot_coefs(train_df)) %>% 
  reduce(bind_rows)

# All coefs distributed approximately normally.
map2(boot_coefs, names(boot_coefs),
     ~ print(hist(.x, breaks = 30, main = .y)))

make_interval <- function(v) {
  m <- mean(v)
  sd <- sd(v)
  upper <- round(m + 1.96 * sd, 3)
  lower <- round(m - 1.96 * sd, 3)
  paste0("[ ", lower, " - ", upper, " ]")
}

means <- 
  boot_coefs %>% 
  summarize_all(mean) %>% 
  t() %>%
  as_tibble(rownames = "term") %>% 
  rename(boot_mean = V1)

ses <-
  boot_coefs %>% 
  summarize_all(~ 1.96 * sd(.x)) %>% 
  t() %>%
  as_tibble(rownames = "term") %>% 
  rename(boot_se = V1)

intervals <- 
  boot_coefs %>% 
  summarize_all(make_interval) %>% 
  t() %>%
  as_tibble(rownames = "term") %>% 
  rename(boot_interval = V1)

# Boot means nearly identical to lm terms.
# Boot SEs consistently a little over double the lm SEs.
inner_join(broom::tidy(train_fit), 
           broom::tidy(test_fit), by = "term") %>% 
  inner_join(means, by = "term") %>% 
  inner_join(ses, by = "term") %>% 
  inner_join(intervals, by = "term") %>% 
  select(-starts_with("statistic"), -starts_with("p.value")) %>% 
  map_if(is.numeric, ~ round(.x, 3)) %>% 
  as_tibble() %>% 
  mutate(
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "last_a1c_value" ~ "Prior HbA1c",
      term == "mean_bg_full_day" ~ "Mean BG",
      term == "black" ~ "Black",
      term == "sd_bg_full_day" ~ "SD BG",
      TRUE ~ NA_character_
    )
  ) %>% 
  rename(
    Coefficient = term,
    `Tr est` = estimate.x,
    `Tr SE` = std.error.x,
    `Tst est` = estimate.y,
    `Tst SE` = std.error.y,
    `Bt est` = boot_mean,
    `Bt SE` = boot_se,
    `Bt CI` = boot_interval
  ) %>% 
  write_csv("lib/analysis/inference-results-boot.csv")

plot(train_fit)

#####
# d #
#####

full_formula <- 
  a1c_value ~ 1 + mean_bg_full_day + sd_bg_full_day + cv_bg_full_day + 
  percent_very_low_full_day + percent_low_full_day + 
  percent_in_target_range_full_day + 
  percent_in_conservative_target_range_full_day + 
  percent_high_full_day + percent_very_high_full_day + 
  in_target_range_mean_bg_full_day + high_mean_bg_full_day + 
  in_target_range_sd_bg_full_day + high_sd_bg_full_day + 
  in_target_range_cv_bg_full_day + high_cv_bg_full_day + 
  age + age:has_age + black + male + hispanic + 
  last_a1c_value

full_fit <- lm(full_formula, train_df)

summary(full_fit)
summary(train_fit)

plot(full_fit)

#####
# e #
#####

# High correlation between mean blood glucose statistics.
train_df %>% 
  sample_frac(0.1) %>% 
  select(contains("mean")) %>% 
  GGally::ggpairs()

# High correlation between time in range statistics.
train_df %>% 
  sample_frac(0.1) %>% 
  select(contains("percent")) %>% 
  GGally::ggpairs()

cor(train_df$mean_bg_full_day, train_df$sd_bg_full_day)
cor(test_df$mean_bg_full_day, test_df$sd_bg_full_day)
