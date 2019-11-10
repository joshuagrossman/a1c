# Analyzes cgm, a1c, and measurement data.

library(glmnet)
library(caret)
library(tidyverse)
library(RANN)
library(elasticnet)
library(modelr)
library(randomForest)
library(gbm)
library(furrr)
library(fs)
library(rlang)
library(gridExtra)

theme_set(theme_bw())

set.seed(1)

# Each row is 1 a1c measurement across the 5 datasets. Potentially multiple a1c
# measurements for each patient.
data <- read_csv("data/features.csv")

# Data contains 5,678 a1c measurements.
nrow(data)

# Many a1c measurements appear to have less than a week of associated cgm data.
# Studies desire at least a week of CGM data as a lead-in or follow-up to an a1c
# measurement, typically.
table(data$cgm_days)

# Sparse data below rounded a1c of 6 mg/dL and above 11 mg/dL. 
# >5.7% a1c indicates prediabetes. >6.5% indicates diabetes.
table(round(data$a1c_value))

# Arbitrarily choosing 5 days as minimum amount of CGM data required allows
# us to capture additional ~1300 a1c measurements.
# Use 5.5% as a minimum and 11.5% as a maximum to allow clean bucketing into
# rounded a1cs.
feature_df <- data %>% 
  filter(! is.na(cgm_days), 
         cgm_days >= 5, 
         between(a1c_value, 5.5, 11.5)) 

# Reduced to 4,415 a1c measurements.
nrow(feature_df)
  
# In Protocol_H study, there are sometimes multiple HbA1c recordings on the same
# day. This was likely to calibrate the a1c measurement. Taking a random one
# from each day where there are multiple recordings.
feature_df <- feature_df %>% 
  group_by(id, a1c_date) %>% 
  sample_n(1) %>% 
  ungroup()

# Reduced to 3,780 a1c measurements.
nrow(feature_df)

# CMetformin and FLEX studies were conducted on individuals with poor diabetes
# control, so it makes sense that mean BGs are higher for those studies.
# Protocol_F falls somewhere in the middle.
feature_df %>% 
  group_by(id, dataset) %>% 
  summarize(mean_bg = mean(mean_bg_full_day, na.rm = T)) %>% 
  ggplot(aes(x = mean_bg, color = dataset)) +
  geom_density()

# CMetformin and FLEX studies were conducted on individuals with poor diabetes
# control, so it makes sense that mean SD BGs are higher for those studies.
# Again, Protocol_F falls somewhere in the middle.
feature_df %>% 
  group_by(id, dataset) %>% 
  summarize(mean_sd_bg = mean(sd_bg_full_day, na.rm = T)) %>% 
  ggplot(aes(x = mean_sd_bg, color = dataset)) +
  geom_density()

# CMetformin is the smallest study by far. Good to keep in mind for 
# cross-dataset validation later on.
feature_df %>% 
  group_by(dataset) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = dataset, y = n)) +
  geom_col() + 
  coord_flip()
  
  mutate(race = ifelse(race == "Black/African American",
                       "Black",
                       race),
         ethnicity = ifelse(ethnicity == "Non-Hispanic",
                            "Not Hispanic or Latino",
                            ethnicity)
  ) %>% 
  ungroup()

