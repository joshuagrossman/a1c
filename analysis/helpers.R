"Helper functions for analysis of CGM data."

library(tidyverse)

load_and_filter <- function(data, write = FALSE) {
  # Loads raw feature data and filters/cleans.
  
  # Don't need some columns.
  # CGM name is usually the same for each study, so study name captures that info.
  # Daytime and nighttime measurements are overly granular.
  feature_df <- select(data, -c("first_day_cgm", "last_day_cgm", "datafile", 
                                "cgm_name", "study_start_date", 
                                "years_since_diag", "metformin"),
                       -ends_with("daytime"),
                       -ends_with("nighttime"))
  
  # A1c of 7-7.5 is a safe target for diabetics.
  feature_df <- mutate(feature_df, 
                       healthy_a1c = as.integer(a1c_value < 7.5))
  
  # # Data contains 5,678 a1c measurements.
  # nrow(data)
  
  # # Many a1c measurements appear to have less than a week of associated cgm data.
  # # Studies desire at least a week of CGM data as a lead-in or follow-up to an a1c
  # # measurement, typically.
  # table(data$cgm_days)
  # 
  # # Sparse data below rounded a1c of 6 mg/dL and above 11 mg/dL. 
  # # >5.7% a1c indicates prediabetes. >6.5% indicates diabetes.
  # table(round(data$a1c_value))
  
  # Arbitrarily choosing 5 days as minimum amount of CGM data required allows
  # us to capture additional ~1300 a1c measurements.
  # Use 5.5% as a minimum and 11.5% as a maximum to allow clean bucketing into
  # rounded a1cs.
  feature_df <- feature_df %>% 
    filter(! is.na(cgm_days), 
           cgm_days >= 5, 
           between(a1c_value, 5.5, 11.5)) 
  
  # # Reduced to 4,415 a1c measurements.
  # nrow(feature_df)
  
  # In Protocol_H study, there are sometimes multiple HbA1c recordings on the same
  # day. This was likely to calibrate the a1c measurement. Taking a random one
  # from each day where there are multiple recordings.
  feature_df <- feature_df %>% 
    group_by(id, a1c_date) %>% 
    sample_n(1) %>% 
    ungroup()
  
  # # Reduced to 3,780 a1c measurements.
  # nrow(feature_df)
  
  # # Some missing data in the CGM statistics. < 2%. Imputation from other CGM
  # # columns seems reasonable. 
  # # 20% of ages are missing. Will use interaction w/ missing indicator.
  # map_dbl(feature_df, ~ mean(is.na(.)))
  
  feature_df <- feature_df %>% 
    mutate(has_age = if_else(! is.na(age), TRUE, FALSE),
           # this looks like mean imputation, but just a placeholder
           age = if_else(has_age, age, mean(age, na.rm = T)))
  
  # # Gender binarized.
  # table(feature_df$gender, useNA = "always")
  # 
  # # Some ages have been imputed beforehand and are non-integers. 
  # table(round(feature_df$age, 1))
  # 
  # # Races and ethnicity not consistently coded.
  # table(feature_df$race)
  # table(feature_df$ethnicity)
  
  # Not enough individuals in races other than White/Black.
  feature_df <- feature_df %>% 
    mutate(race = if_else(race == "Black/African American", "Black", race),
           ethnicity = if_else(ethnicity == "Non-Hispanic", 
                               "Not Hispanic or Latino",
                               ethnicity)) %>% 
    filter(race %in% c("Black", "White"),
           ethnicity %in% c("Hispanic or Latino", "Not Hispanic or Latino"))
  
  # Recoding for ease later on
  feature_df <- feature_df %>% 
    mutate(black = if_else(race == "Black", 1, 0),
           male = if_else(gender == "M", 1, 0),
           hispanic = if_else(ethnicity == "Hispanic or Latino", 1, 0)) %>% 
    select(-race, -gender, -ethnicity)
  
  # # RT_CGM used datetimes instead of date, and supplied a weird origin.
  # # Don't need exact dates, just care about differences between dates, so
  # # using default origin since their origin cannot be determined.
  feature_df <- feature_df %>%
    mutate(is_datetime = a1c_date > 30000,
           a1c_date = if_else(is_datetime,
                              as_date(as_datetime(a1c_date)),
                              as_date(a1c_date))) %>%
    select(-is_datetime)
    
  
  # # Reduced to 3,559 a1c measurements.
  # nrow(feature_df)
  
  if (write) {
    write_csv(feature_df, "data/cleaned_features.csv")
  }
  
  feature_df
}

extract_last_a1c <- function(feature_df) {
  grouped_df <- feature_df %>%
    group_by(id) %>% 
    arrange(a1c_date, .by_group = T)
  
  new_cols <- grouped_df %>% 
    group_map(~ select(., last_a1c_value = a1c_value, 
                       last_a1c_date = a1c_date,
                       last_a1c_resid = a1c_resid,
                       last_healthy_a1c = healthy_a1c,
                       last_healthy_resid = healthy_resid)) %>% 
    map(~ head(., -1)) %>% 
    map(~ bind_rows(tibble(last_a1c_value = NA, 
                           last_a1c_date = NA,
                           last_a1c_resid = NA,
                           last_healthy_a1c = NA,
                           last_healthy_resid = NA), .)) %>% 
    reduce(bind_rows)
  
  bind_cols(grouped_df, new_cols) %>% 
    mutate(days_since_last_a1c = a1c_date - last_a1c_date) %>% 
    select(-last_a1c_date) %>% 
    ungroup()
}
  