# Constants
# DONE Max number of records for each day
# DONE Low BG limit
# DONE High BG limit
# DONE ranges, sleep 12 AM - 6 AM, wake 6 AM - 12 AM, 24 hours
# DONE 70-80% data required (make sure satisfied if 15 min recordings)


# Minimum BG cutoff (check data?)
# Maximum BG cutoff (check data?)


# Features
# DONE Count of valid records at each day interval
# DONE Fraction of max records valid at each day interval
# DONE Average BG at each day interval
# DONE TIR in each day interval 70-180, 70-140
# DONE SD at each interval
# DONE Time very high >250
# DONE Time high >180 <250
# DONE Time low >54 <70
# DONE Time very low <54
# DONE Area under the curve - average at each day interval, normalized by day

# Number of hypo / hyperglycemic events per day
# CV at each interval

# DONE How to filter missing data? By day interval? Or by day? -- Full Day in the past

# DONE How to average/sd results? Over all measurements (unweighted) or by day? -- Average by day

# DONE How to calculate AUC for missing data? Weight it by the data points you
# actually have

# DONE How to report AUC feature? Average AUC?  -- Normalize it by day

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

source("/lib/load.R")

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

MIN_DATA_REQUIRED <- 0.7

################################## CALCULATE CGM FEATURES ######################

make_cgm_feature_df <- function(cgm_data_with_id) {
  # Input: List containing CGM data for a particular patient and the patient's id
  #
  # Output: One-row dataframe containing summary features for the patient's CGM data.
  
  id <- as.character(cgm_data_with_id$id)
  cgm_data <- cgm_data_with_id$data
  
  sufficient_cgm_data <- filter_insufficient_data(cgm_data)
  
  calculate_all_bg_stats(sufficient_cgm_data)
}

############# FILTER TOO OLD, TOO RECENT, AND INSUFFICIENT DATA ################

Mode <- function(x) {
  
  # R doesn't have a base mode function...
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

calculate_minutes_per_record <- function(bg_df, n = 1000) {
  # Estimates the amount of time between CGM readings for a given patient's
  # CGM data. 
  #
  # Assumes that all CGM data for a particular patient will have one
  # common time between recordings. 
  
  datetimes <- bg_df$datetime
  
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
  
make_intervals <- function(bg_df) {
  # Extracts date and hour from timestamps, and indicates 
  # measurements that occur during nighttime. 
  
  bg_df %>% 
    mutate(cgm_day = lubridate::date(datetime),
           cgm_hour = lubridate::hour(datetime),
           nighttime = cgm_hour < DAY_START)
}

calculate_records_per_hour <- function(minutes_per_record) {
  # Determine max possible CGM recordings per hour.
  
  if (! minutes_per_record %in% MINUTES_BETWEEN_RECORDS) {
    print(minutes_per_record)
    abort("CGM data not recorded in 5, 10, 15, or 20 minute intervals.")
  }
  
  60 / minutes_per_record
}

filter_data_by_date <- function(bg_df, 
                                start_date_chr = "1900-01-01",
                                end_date_chr = "2100-12-31") {
  # Removed CGM data for invalid dates.
  
  start_date <- as.Date(start_date_chr)
  end_date <- as.Date(end_date_chr)
  
  bg_df <-
    filter(bg_df, between(datetime, start_date, end_date))
}
                                
filter_insufficient_data <- function(bg_df) {
  # Counts total possible CGM records for a given day and calculates the fraction
  # of records actually recorded. Removes data from days with less than 
  # MIN_DATA_REQUIRED fraction of records present.
  
  cgm_records_per_hour <- calculate_minutes_per_record(bg_df) %>% 
                            calculate_records_per_hour()
  
  interval_df <- make_intervals(bg_df)
  
  exclude_df <-
    interval_df %>%
    group_by(cgm_day) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(prop = n / (cgm_records_per_hour * 24)) %>%
    filter(prop < MIN_DATA_REQUIRED) %>%
    select(cgm_day)
  
  anti_join(interval_df, exclude_df, by = c("cgm_day"))
}

################################## SIMPLE BG STATS #############################

calculate_cgm_days <- function(bg_df) {
  # Calculates the total number of days for which CGM data is available
  
  c(cgm_days = length(unique(bg_df$cgm_day)))
}

calculate_simple_bg_stat <- function(bg_df, f, stat_name) {
  # bg_df: dataframe with CGM data, 
  # f: function to perform on glucose values,
  # stat_name: a string describing f
  #
  # Output: a named vector of the values calculated by applying f to glucose
  # values during the full day, nighttime, and daytime.
  
  bg <- bg_df$glucose
  nighttime <- bg_df$nighttime
  bg_nighttime <- bg[nighttime]
  bg_daytime <- bg[! nighttime]
  
  bg_stat <- c(f(bg), f(bg_nighttime), f(bg_daytime))
  
  names(bg_stat) <- c(str_c("full_day_", stat_name, "_bg"),
                      str_c("nighttime_", stat_name, "_bg"),
                      str_c("daytime_", stat_name, "_bg"))
  
  bg_stat
}

################################## TIME IN RANGE ###############################

make_bg_booleans <- function(bg_df) {
  # Determines which range a particular CGM recording falls in.
  
  bg_df %>% 
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

calculate_percent_in_given_range <- function(bg_df, range_name) {
  # Calculates the percent of valid CGM readings that fall in the supplied
  # range_name (e.g. "very_high")
  #
  # Assumes that range_name is also a column in bg_df (created by make_bg_booleans)
  
  if (! range_name %in% colnames(bg_df)) {
    abort(str_c(range_name, " column not contained in argument `bg_df` of `calculate_percent_in_given_range`."))
  }
  
  in_range <- pull(bg_df, range_name)
  full_day_percent_in_range <- na_mean(in_range)
  
  nighttime <- bg_df$nighttime
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

################################## AUC #########################################

calculate_auc <- function(bg_df) {
  # Calculates the AUC for a given day of CGM data, and scales the AUC by the
  # number of valid recordings for that day
  
  time_diffs <- bg_df$datetime %>% 
                  diff() %>% 
                  time_length(unit = "minutes") %>% 
                  abs() %>% 
                  round()
  
  minutes_per_record <- calculate_minutes_per_record(bg_df)
  
  max_records_per_day <- 24 * calculate_records_per_hour(minutes_per_record)
  
  glucose <- bg_df$glucose
  
  # if time_diff is the same as normal recording interval, calculate auc, else 0
  aucs <- pmap(list(head(glucose, -1),
                    glucose[-1],
                    time_diffs),
       ~ if (..3 == minutes_per_record) {..3 * (..1 + ..2) / 2} else 0)
  
  aucs <- unlist(aucs)
  
  n_valid_aucs <- sum(aucs != 0)
  
  # if 288 possible recordings in a day, 287 possible aucs can be recorded
  sum(aucs) * ((max_records_per_day - 1)  / n_valid_aucs)
}

calculate_average_auc <- function(bg_df, range_name = NULL) {
  
  # if (! is.null(range_name)) {
  #   
  #   if (! range_name %in% colnames(bg_df)) {
  #     abort(str_c(range_name, " column not contained in argument `bg_df` of `calculate_average_auc`."))
  #   }
  #   
  #   in_range <- pull(bg_df, range_name)
  #   
  #   minutes_per_record <- calculate_minutes_per_record(bg_df)
  #   
  #   records_per_half_hour <- ceiling(60 / minutes_per_record)
  #   
  #   # collapse in_range 0 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 into 0 2 000 3 000 4
  #   # if the numbers are >= records_per_half_hour, expand into corresponding number of 1s
  #   # else expand into corresponding numbers of 0s
  #
  #   in_range_at_least_half_hour <- 
  #   
  #   bg_df <-
  #     filter(bg_df, in_range_at_least_half_hour)
  # }
  
  average_auc <- 
    bg_df %>% 
    group_by(cgm_day) %>% 
    group_map(~ calculate_auc(.x)) %>% 
    unlist() %>% 
    mean()
  
  names(average_auc) <- str_c("auc_", range_name)
  
  average_auc
}

calculate_all_bg_stats <- function(bg_df) {
  
  bg_df <- make_bg_booleans(bg_df)
  
  features <- c(calculate_cgm_days(bg_df),
                calculate_simple_bg_stat(bg_df, na_mean, "mean"),
                calculate_simple_bg_stat(bg_df, na_sd, "sd"),
                calculate_average_auc(bg_df))
  
  in_range_stats <- map(RANGE_NAMES, 
                        ~ calculate_percent_in_given_range(bg_df, .))
  
  # auc_stats <- map(RANGE_NAMES, 
  #                       ~ calculate_average_auc(bg_df, .))
  
  features <- c(features, unlist(in_range_stats))#, unlist(auc_stats))
}



################################## DEPRECATED ##################################

# load_individual_file <- function(filepath) {
#   
#   if (path_ext(filepath) != "csv") {
#     abort("CGM datafile must be a CSV.")
#   }
#   
#   bg_df <- read_csv(filepath)
#   
#   if (! all(colnames(bg_df) %in% CGM_HEADERS)) {
#     abort(str_c("CGM datafile must contain ", CGM_HEADERS))
#   }
#   
#   id <- path_ext_remove(path_file(filepath))
#   
#   list(id = id, data = bg_df)
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