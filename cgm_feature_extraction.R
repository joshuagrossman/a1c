# Constants
# Max number of records for each day
# Low BG limit
# High BG limit
# Minimum BG cutoff (check data)
# Maximum BG cutoff (check data)
# ranges, sleep 12 AM - 6 AM, wake 6 AM - 12 AM, 24 hours
# 70-80% data required (make sure satisfied if 15 min recordings)


# Features
# Count of valid records at each day interval
# Fraction of max records valid at each day interval
# Average BG at each day interval
# TIR in each day interval 70-180, 70-140
# SD at each interval
# CV at each interval
# Time very high >250
# Time high >180 <250
# Time low >54 <70
# Time very low <54
# Area under the curve
# Number of hypo / hyperglycemic events per day

################################## LIBRARIES ###################################

library(tidyverse)
library(fs)
library(rlang)
library(lubridate)

################################## GLOBALS #####################################

MINUTES_BETWEEN_RECORDS <- c(5, 15)

DAY_START <- 6

LOW_BG_LIMIT <- 70
VERY_LOW_BG_LIMIT <- 54
HIGH_BG_LIMIT <- 180
VERY_HIGH_BG_LIMIT <- 250

MIN_BG_CUTOFF <- 40
MAX_BG_CUTOFF <- 400

MIN_DATA_REQUIRED <- 0.7

CGM_HEADERS <- c("glucose", "datetime")

################################## LOAD DATA ###################################

load_file <- function(filepath) {
  
  if (path_ext(filepath) != "csv") {
    abort("CGM datafile must be a CSV.")
  }
  
  df <- read_csv(filepath)
  
  if (! all(colnames(df) %in% CGM_HEADERS)) {
    abort(str_c("CGM datafile must contain ", CGM_HEADERS))
  }
  
  id <- path_ext_remove(path_file(filepath))
  
  list(id = id, data = df)
}

load_dir <- function(dirpath) {
  
  if (! is_dir(dirpath)) {
    abort("CGM data folder must be a directory.")
  }
  
  dir_map(dirpath, load_file)
}


################################## FILTER INSUFFICIENT DATA ####################

calculate_minutes_per_record <- function(df) {
  df$datetime %>% 
    head(1000) %>% 
    diff() %>% 
    time_length(unit = "minutes") %>% 
    median(na.rm = T)
}
  
make_intervals <- function(df) {
  df %>% 
    mutate(cgm_day = lubridate::date(datetime),
           cgm_hour = lubridate::hour(datetime),
           nighttime = cgm_hour < DAY_START)
}

count_records_per_interval <- function(df) {
  
  if (! calculate_minutes_per_record(df) %in% MINUTES_BETWEEN_RECORDS) {
    abort("CGM data not recorded in 5 or 15 minute intervals.")
  }
  
  cgm_record_multiplier <- 1
  
  if (median_time_between_records == 5) {
    cgm_record_multiplier <- 3
  }
  
  cgm_records_per_hour <- 4 * cgm_record_multiplier
    
  daytime_record_multiplier <- ((24 - DAY_START) / DAY_START) - 1
  
  interval_df <- make_intervals(df)
    
  exclude_df <- 
    interval_df %>% 
    group_by(cgm_day, nighttime) %>%
    summarize(n = n()) %>% 
    ungroup() %>% 
    mutate(prop = n / (cgm_records_per_hour * 
                            (DAY_START * 
                               (1 + (! nighttime) * 
                               daytime_record_multiplier)))) %>% 
    filter(prop < MIN_DATA_REQUIRED) %>% 
    select(cgm_day, nighttime)
  
  anti_join(interval_df, exclude_df, by = c("cgm_day", "nighttime"))
}

################################## BOOLEANS FOR BG STATS #######################

make_bg_booleans <- function(df) {
  
  bg_df <-
    df %>% 
    mutate(very_low <- glucose < VERY_LOW_BG_LIMIT,
           low <- glucose < LOW_BG_LIMIT,
           in_target_range <- between(glucose, LOW_BG_LIMIT, HIGH_BG_LIMIT),
           high <- glucose > HIGH_BG_LIMIT,
           very_high <- glucose > VERY_HIGH_BG_LIMIT)
  
  list(data = bg_df, features = NULL)
}

################################## CALCULATE BG STATS #####################################

calculate_mean_bg <- function(d) {
  
  df <- d$data
  features <- d$features
  
  features <- c(features, 
                mean_bg = mean(df$glucose, na.rm = T))
  
  list(data = df, features = features)
}

calculate_percent_in_range <- function(df, range_name) {
  
  in_range <- `$`(df, range_name)
  full_day_percent_in_range <- mean(in_range, na.rm = T)
  
  nighttime <- df$nighttime
  in_range_nighttime <- in_range[nighttime]
  in_range_daytime <- in_range[! nighttime]
    
  nighttime_percent_in_range <- mean(in_range_nighttime, na.rm = T)
  daytime_percent_in_range <- mean(in_range_daytime, na.rm = T)
  
  percent_in_range_features <- 
    c(full_day_percent_in_range,
      nighttime_percent_in_range,
      daytime_percent_in_range)
  
  names(percent_in_range_features) <- 
    c(str_c("percent_", range_name, "_full_day"),
      str_c("percent_", range_name, "_nighttime"),
      str_c("percent_", range_name, "_daytime"))
  
  percent_in_range_features
}

