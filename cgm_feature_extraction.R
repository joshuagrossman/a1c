# Constants
# DONE Max number of records for each day
# DONE Low BG limit
# DONE High BG limit
# Minimum BG cutoff (check data)
# Maximum BG cutoff (check data)
# ranges, sleep 12 AM - 6 AM, wake 6 AM - 12 AM, 24 hours
# DONE 70-80% data required (make sure satisfied if 15 min recordings)


# Features
# DONE Count of valid records at each day interval
# DONE Fraction of max records valid at each day interval
# DONE Average BG at each day interval
# DONE TIR in each day interval 70-180, 70-140
# DONE SD at each interval
# CV at each interval
# DONE Time very high >250
# DONE Time high >180 <250
# DONE Time low >54 <70
# DONE Time very low <54
# Area under the curve - average at each day interval?
# Number of hypo / hyperglycemic events per day

# How to filter missing data? By day interval? Or by day? -- Full Day in the past

# How to average/sd results? Over all measurements (unweighted) or by day? -- Average by day

# How to calculate AUC for missing data? Weight it by the data points you
# actually have

# How to report AUC feature? Average AUC?  -- Normalize it by day

# Brainstorming additional features worth extracting

# AUC at high / low values 

# Rate of change - only include times of 1-2 hours of increase

# Rate of change of AUC across hours, versus half hours, versus 15 minutes... as 
# a way to define the complexity of the curve

# What is the 90th percentile of the highs for a particular interval

# Percentile statistics as opposed to average (don't use max/min)

# Google definition of fractal dimensions for curves

# Keeping a list of how to analyze CGM data different ways - averaging over days versys averaging over all points

# Sensitivity analysis - calculate all of these metrics for a particular person
# removing 1, 2, 3, ... data points at random, quantify bias of missing data

# Random forest, LASSO, GBM, and XGBoost, compare against eA1c

# 10-fold CV to tune parameters with grid search, good metric to report is R^2, 
# scatterplot for both models of actual versus predicted, color the points by the
# value of the demographic of interest

################################## LIBRARIES ###################################

library(tidyverse)
library(fs)
library(rlang)
library(lubridate)
library(furrr)

plan(multisession(workers = availableCores()))

################################## GLOBALS #####################################

MINUTES_BETWEEN_RECORDS <- c(5, 10, 15, 20)

DAY_START <- 6
# DAY_END assumed as midnight. Potentially update in the future.

LOW_BG_LIMIT <- 70
VERY_LOW_BG_LIMIT <- 54
CONSERVATIVE_HIGH_BG_LIMIT <- 140
HIGH_BG_LIMIT <- 180
VERY_HIGH_BG_LIMIT <- 250

RANGE_NAMES <- c("very_low", 
                 "low", 
                 "in_target_range", 
                 "in_conservative_target_range", 
                 "high", 
                 "very_high")

MIN_BG_CUTOFF <- 40
MAX_BG_CUTOFF <- 400

MIN_DATA_REQUIRED <- 0.7

CGM_HEADERS <- c("glucose", "datetime")

################################## LOAD DATA ###################################

load_large_file <- function(filepath) {
  # Input: Filepath to a CSV containing CGM data for more than one patient. The
  # CSV must contain patient id info, glucose values, and timestamps.
  #
  # Output: A list of lists, with each list containing the patient_id and all
  # associated CGM data.
  
  if (path_ext(filepath) != "csv") {
    abort("CGM datafile must be a CSV.")
  }
  
  df <- read_csv(filepath)
  
  headers_with_id <- c("id", CGM_HEADERS)
  
  if (! all(colnames(df) %in% headers_with_id)) {
    abort(str_c("CGM datafile must contain ", headers_with_id))
  }
  
  group_by(df, id) %>% 
    group_map(~ list(id = .y, data = .x))
}

# load_individual_file <- function(filepath) {
#   
#   if (path_ext(filepath) != "csv") {
#     abort("CGM datafile must be a CSV.")
#   }
#   
#   df <- read_csv(filepath)
#   
#   if (! all(colnames(df) %in% CGM_HEADERS)) {
#     abort(str_c("CGM datafile must contain ", CGM_HEADERS))
#   }
#   
#   id <- path_ext_remove(path_file(filepath))
#   
#   list(id = id, data = df)
# }
# 
# load_dir <- function(dirpath) {
#   
#   if (! is_dir(dirpath)) {
#     abort("CGM data folder must be a directory.")
#   }
#   
#   dir_map(dirpath, load_individual_file)
# }

################################## PIPELINE ####################################

make_feature_df <- function(filepath) {
  # Input: Filepath to a CSV containing CGM data for all patients.
  #
  # Output: Dataframe containing summary features for each patient's CGM data. 
  
  patient_data_with_ids <- load_large_file(filepath)
  
  ids <- map_chr(patient_data_with_ids, ~ as.character(.x$id))
  patient_data <- map(patient_data_with_ids, ~ .x$data)
  
  augmented_patient_data <- future_map(patient_data, count_records_per_interval)
  
  future_map(augmented_patient_data, calculate_all_bg_stats) %>% 
    reduce(rbind) %>%  
    as_tibble() %>% 
    mutate(id = ids)
}

################################## FILTER INSUFFICIENT DATA ####################

Mode <- function(x) {
  # R doesn't have a base mode function...
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

calculate_minutes_per_record <- function(df, n = 1000) {
  # Estimates the amount of time between CGM readings for a given patient's
  # CGM data. Assumes that all CGM data for a particular patient will have one
  # common time between recordings. 
  
  datetimes <- df$datetime
  
  n_rows <- length(datetimes)
  
  # select n sequential rows at a random start point
  if (n_rows > n) {
     start <- sample(1:(n_rows - n + 1), 1)
     datetimes <- datetimes[start:(start + n - 1)]
  }
  
  datetimes %>%  
    diff() %>% 
    time_length(unit = "minutes") %>% 
    abs() %>% 
    round() %>% 
    Mode()
}
  
make_intervals <- function(df) {
  # Converts raw timestamps to date and hour equivalents, and indicates 
  # measurements that occur during nighttime. 
  
  df %>% 
    mutate(cgm_day = lubridate::date(datetime),
           cgm_hour = lubridate::hour(datetime),
           nighttime = cgm_hour < DAY_START)
}

calculate_records_per_hour <- function(df) {
  # Determine max possible CGM recordings per hour.
  
  mode_time_between_records <- calculate_minutes_per_record(df)
  
  if (! mode_time_between_records %in% MINUTES_BETWEEN_RECORDS) {
    print(df)
    abort("CGM data not recorded in 5, 10, 15, or 20 minute intervals.")
  }
  
  # calculates how many times more you are recording per hour compared to 
  # 20 minute recording intervals (3 times per hour)
  
  cgm_record_multiplier <- 1
  
  if (mode_time_between_records == 5) {
    cgm_record_multiplier <- 4
  }
  
  if (mode_time_between_records == 10) {
    cgm_record_multiplier <- 2
  }
  
  if (mode_time_between_records == 15) {
    cgm_record_multiplier <- 4/3
  }
  
  3 * cgm_record_multiplier
}

count_records_per_interval <- function(df) {
  # Counts total possible CGM records for a given day and calculates the fraction
  # of records actually recorded. Removes data from days with less than 
  # MIN_DATA_REQUIRED fraction of records present.
  
  cgm_records_per_hour <- calculate_records_per_hour(df)
  
  interval_df <- make_intervals(df)
  
  exclude_df <-
    interval_df %>%
    group_by(cgm_day) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(prop = n / (cgm_records_per_hour * 24)) %>%
    filter(prop < MIN_DATA_REQUIRED) %>%
    select(cgm_day)
  
  anti_join(interval_df, exclude_df, by = c("cgm_day"))
    
  # DEPRECATED: this code removed data by day/night interval
  #.
  # TODO: Make this more robust with DAY_START and DAY_END instead of using
  # midnight default as day end
  # 
  # daytime_record_multiplier <- ((24 - DAY_START) / DAY_START) - 1
  #
  # exclude_df <-
  #   interval_df %>%
  #   group_by(cgm_day, nighttime) %>%
  #   summarize(n = n()) %>%
  #   ungroup() %>%
  #   mutate(prop = n / (cgm_records_per_hour *
  #                       (DAY_START *
  #                         (1 + (! nighttime) * daytime_record_multiplier)))) %>%
  #   filter(prop < MIN_DATA_REQUIRED) %>%
  #   select(cgm_day, nighttime)
  # 
  # anti_join(interval_df, exclude_df, by = c("cgm_day", "nighttime"))
}

################################## CALCULATE BG STATS ##########################

calculate_cgm_days <- function(df) {
  # Calculates the total number of days for which CGM data is available
  
  c(cgm_days = length(unique(df$cgm_day)))
}

calculate_simple_bg_stat <- function(df, f, stat_name) {
  # Input: df: dataframe with CGM data, f: function to perform on glucose values,
  # stat_name: a string describing f
  #
  # Output: a named vector of the values calculated by applying f during the
  # full day, nighttime, and daytime.
  
  bg <- df$glucose
  nighttime <- df$nighttime
  bg_nighttime <- bg[nighttime]
  bg_daytime <- bg[! nighttime]
  
  bg_stat <- c(f(bg), f(bg_nighttime), f(bg_daytime))
  
  names(bg_stat) <- c(str_c("full_day_", stat_name, "_bg"),
                      str_c("nighttime_", stat_name, "_bg"),
                      str_c("daytime_", stat_name, "_bg"))
  
  bg_stat
}

make_bg_booleans <- function(df) {
  # Determines which range a particular CGM recording falls in.
  
  df %>% 
    mutate(very_low = glucose < VERY_LOW_BG_LIMIT,
           low = glucose < LOW_BG_LIMIT,
           in_target_range = between(glucose, LOW_BG_LIMIT, HIGH_BG_LIMIT),
           in_conservative_target_range = between(glucose, LOW_BG_LIMIT, CONSERVATIVE_HIGH_BG_LIMIT),
           high = glucose > HIGH_BG_LIMIT,
           very_high = glucose > VERY_HIGH_BG_LIMIT)
}

na_mean <- function(v) {
  mean(v, na.rm = T)
}

na_sd <- function(v) {
  sd(v, na.rm = T)
}

calculate_percent_in_given_range <- function(df, range_name) {
  # Calculates the percent of valid CGM readings that fall in the supplied
  # range_name (e.g. "very_high")
  #
  # Assumes that range_name is also a column in df (created by make_bg_booleans)
  
  in_range <- pull(df, range_name)
  full_day_percent_in_range <- na_mean(in_range)
  
  nighttime <- df$nighttime
  in_range_nighttime <- in_range[nighttime]
  in_range_daytime <- in_range[! nighttime]
    
  nighttime_percent_in_range <- na_mean(in_range_nighttime)
  daytime_percent_in_range <- na_mean(in_range_daytime)
  
  percent_in_range_features <- 
    c(full_day_percent_in_range,
      nighttime_percent_in_range,
      daytime_percent_in_range)
  
  names(percent_in_range_features) <- 
    c(str_c(range_name, "_full_day"),
      str_c(range_name, "_nighttime"),
      str_c(range_name, "_daytime"))
  
  percent_in_range_features
}

calculate_auc <- function(df) {
  # Calculates the AUC for a given day of CGM data, and scales the AUC by the
  # number of valid recordings for that day
  
  time_diffs <- df$datetime %>% 
                  diff() %>% 
                  time_length(unit = "minutes") %>% 
                  abs() %>% 
                  round()
  
  n_records <- calculate_cgm_days(df)
  
  glucose <- df$glucose
  
  # TODO: if time_diff is the same as normal recording interval, calculate, else ignore
  aucs <- pmap(list(head(glucose, -1),
                    glucose[-1],
                    time_diffs),
       ~ ..3 * (..1 + ..2) / 2)
  
  sum(unlist(aucs)) / n_records
}

calculate_average_auc <- function(df) {
  average_auc <- 
    df %>% 
    group_by(cgm_day) %>% 
    group_map(~ calculate_auc(.x)) %>% 
    unlist() %>% 
    mean()
    
  c(average_auc = average_auc)
}

calculate_all_bg_stats <- function(df) {
  
  bg_df <- make_bg_booleans(df)
  
  features <- c(calculate_cgm_days(bg_df),
                calculate_simple_bg_stat(bg_df, na_mean, "mean"),
                calculate_simple_bg_stat(bg_df, na_sd, "sd"),
                calculate_average_auc(bg_df))
  
  in_range_stats <- map(RANGE_NAMES, 
                        ~ calculate_percent_in_given_range(bg_df, .))
  
  features <- c(features, unlist(in_range_stats))
}