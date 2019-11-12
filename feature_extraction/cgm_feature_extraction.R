# Extracts features from CGM data.

source("lib/feature_extraction/load.R")

################################## GLOBALS #####################################

MINUTES_BETWEEN_RECORDS <- c(5, 10, 15, 20)

# Standards based on Danne et al., 2017 (International Consensus on Use of
# Continuous Glucose Monitoring)

DAY_START <- 6
# DAY_END assumed as midnight.

# TODO: Figure out if CGM times are recorded as local time or GMT

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

# most relevant sources of variation in BG data
CORE_RANGE_NAMES <- c("in_target_range", 
                      "high")

# minimum fraction of possible readings
MIN_DATA_REQUIRED <- 0.7

# data points required to calculate summary stats
MIN_RECORDINGS_REQUIRED <- 10

################################## CALCULATE CGM FEATURES ######################

make_cgm_feature_df <- function(cgm_data_with_id,
                                data_name = "data") {
  # Input: List containing CGM $data for a particular patient and the patient's
  # $id, as generate by `load_file()`
  #
  # Output: One-row dataframe containing summary features for the patient's CGM
  # data.
  
  cgm_data <- cgm_data_with_id[[data_name]]
  
  if (! is_tibble(cgm_data)) {
    warning("No cgm data provided. Returning NULL.")
    return(NULL)
  }
  
  if (nrow(cgm_data) < MIN_RECORDINGS_REQUIRED) {
    warning(str_c("Less than ", MIN_RECORDINGS_REQUIRED, 
                  " data points. Returning NULL."))
    return(NULL)
  }
  
  cgm_data <- arrange(cgm_data, datetime)
  
  sufficient_cgm_data <- filter_insufficient_data(cgm_data)
  
  if (is.null(sufficient_cgm_data)) {
    # warning message is supplied by `filter_insufficient_data` function
    return(NULL)
  }
  
  calculate_all_bg_stats(sufficient_cgm_data)
}

################## CGM FEATURE EXTRACTION PIPELINE #############################

calculate_all_bg_stats <- function(bg_df) {
  # Applies all bg stats functions to bg_df
  
  if (! is_tibble(bg_df)) {
    abort("Argument `bg_df` of `calculate_all_bg_stats` is not a data frame.")
  }
  
  bg_df <- make_bg_booleans(bg_df)
  
  features <- c(calculate_cgm_days(bg_df),
                determine_first_and_last_cgm_day(bg_df),
                calculate_simple_bg_stat(bg_df, na_mean, "mean"),
                calculate_simple_bg_stat(bg_df, na_sd, "sd"),
                calculate_simple_bg_stat(bg_df, na_cv, "cv"))
  
  percent_in_range <- map(RANGE_NAMES, 
                          ~ calculate_percent_in_given_range(bg_df, .))
  
  # only need summary stats at the main ranges, not enough data at very low etc.
  mean_in_range <- map(CORE_RANGE_NAMES, 
                       ~ calculate_stat_in_given_range(bg_df, 
                                                       ., 
                                                       na_mean, 
                                                       "mean"))
  
  sd_in_range <- map(CORE_RANGE_NAMES, 
                     ~ calculate_stat_in_given_range(bg_df, 
                                                     ., 
                                                     na_sd, 
                                                     "sd"))
  
  cv_in_range <- map(CORE_RANGE_NAMES, 
                     ~ calculate_stat_in_given_range(bg_df, 
                                                     ., 
                                                     na_cv, 
                                                     "cv"))
  
  c(features, unlist(c(percent_in_range, 
                       mean_in_range, 
                       sd_in_range,
                       cv_in_range)))
}


############# FILTER TOO OLD, TOO RECENT, AND INSUFFICIENT DATA ################

Mode <- function(x) {
  # R doesn't have a base mode function...
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

calculate_minutes_per_record <- function(bg_df, 
                                         n = 1000,
                                         datetime_name = "datetime") {
  # Estimates the amount of time between CGM readings for a given patient's
  # CGM data. 
  #
  # Assumes that all CGM data for a particular patient will have one
  # common time between all recordings (e.g., patient didn't switch CGMs)
  
  datetimes <- pull(bg_df, datetime_name)
  
  if (is_character(datetimes)) {
    abort("Datetimes parsed as strings. Convert to datetime and re-run.")
  }
  
  n_rows <- length(datetimes)
  
  # don't use NA datetimes to calculate time between records
  datetimes <- datetimes[! is.na(datetimes)]
  
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
  
make_intervals <- function(bg_df, 
                           datetime_name = "datetime") {
  # Extracts date and hour from timestamps, and indicates 
  # measurements that occur during nighttime. 
  
  bg_df %>% 
    mutate(cgm_day = lubridate::date(.data[[datetime_name]]),
           cgm_hour = lubridate::hour(.data[[datetime_name]]),
           nighttime = cgm_hour < DAY_START)
}

filter_data_by_date <- function(bg_df, 
                                start_date_chr = "1900-01-01",
                                end_date_chr = "2100-12-31",
                                datetime_name = "datetime") {
  # Removes CGM data for invalid dates.
  
  start_date <- as_date(start_date_chr)
  end_date <- as_date(end_date_chr)
  
  bg_df <-
    filter(bg_df, between(as_date(.data[[datetime_name]]), 
                          start_date, 
                          end_date))
}
                                
filter_insufficient_data <- function(bg_df, 
                                     cgm_day_name = "cgm_day") {
  # Counts total possible CGM records for a given day and calculates the fraction
  # of possible records actually recorded. Removes data from days with less than 
  # MIN_DATA_REQUIRED fraction of records present.
  
  minutes_per_record <- calculate_minutes_per_record(bg_df)
  
  if (! minutes_per_record %in% MINUTES_BETWEEN_RECORDS) {
    warning(str_c("CGM data not recorded consistently in 5, 10, 15, or ",
                  "20 minute intervals. Returning NULL."))
    return(NULL)
  }
  
  cgm_records_per_hour <- 60 / minutes_per_record
  
  interval_df <- make_intervals(bg_df)
  
  exclude_df <-
    interval_df %>%
    group_by(cgm_day) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(prop = n / (cgm_records_per_hour * 24)) %>%
    filter(prop < MIN_DATA_REQUIRED) %>%
    select(cgm_day)
  
  anti_join(interval_df, exclude_df, by = c(cgm_day_name))
}

################################## SUMMARY BG STATS ############################

na_mean <- function(v) {
  mean(v, na.rm = T)
}

na_sd <- function(v) {
  sd(v, na.rm = T)
}

na_cv <- function(v) {
  na_sd(v) / na_mean(v)
}

check_if_sufficient_length <- function(v) {
  
  if (length(v) < MIN_RECORDINGS_REQUIRED) {
    return(NA)
  }
  v
}

determine_first_and_last_cgm_day <- function(bg_df, 
                                             cgm_day_name = "cgm_day") {
  # Returns the first and last days of CGM data as a vector
  
  cgm_days <- unique(pull(bg_df, cgm_day_name))
  
  c(first_day_cgm = min(cgm_days, na.rm = T),
    last_day_cgm = max(cgm_days, na.rm = T))
}

calculate_cgm_days <- function(bg_df, 
                               cgm_day_name = "cgm_day") {
  # Calculates the total number of days for which CGM data is available
  
  c(cgm_days = pull(bg_df, cgm_day_name) %>% 
                unique() %>% 
                length())
}

calculate_simple_bg_stat <- function(bg_df, f, stat_name) {
  # bg_df: dataframe with CGM data, 
  # f: function to apply to glucose values,
  # stat_name: a string describing f
  #
  # Output: a named vector of the values calculated by applying f to glucose
  # values during the full day, nighttime, and daytime.
  
  bg <- bg_df$glucose
  nighttime <- bg_df$nighttime
  bg_nighttime <- bg[nighttime]
  bg_daytime <- bg[! nighttime]
  
  bg_stat <- c(f(bg %>% check_if_sufficient_length), 
               f(bg_nighttime %>% check_if_sufficient_length), 
               f(bg_daytime %>% check_if_sufficient_length))
  
  names(bg_stat) <- c(str_c(stat_name, "_bg_full_day"),
                      str_c(stat_name, "_bg_nighttime"),
                      str_c(stat_name, "_bg_daytime"))
  
  bg_stat
}

################################## TIME IN RANGE ###############################

make_bg_booleans <- function(bg_df) {
  # Determines which range a particular CGM recording falls in.
  
  bg_df %>% 
    mutate(very_low = glucose < VERY_LOW_BG_LIMIT,
           low = between(glucose,
                         VERY_LOW_BG_LIMIT,
                         LOW_BG_LIMIT),
           in_target_range = between(glucose, 
                                     LOW_BG_LIMIT, 
                                     HIGH_BG_LIMIT),
           in_conservative_target_range = between(glucose, 
                                                  LOW_BG_LIMIT, 
                                                  CONSERVATIVE_HIGH_BG_LIMIT),
           high = between(glucose,
                          HIGH_BG_LIMIT,
                          VERY_HIGH_BG_LIMIT),
           very_high = glucose > VERY_HIGH_BG_LIMIT)
}

calculate_percent_in_given_range <- function(bg_df, range_name) {
  # Calculates the percent of valid CGM readings that fall in the supplied
  # range_name (e.g. "very_high") for each day interval
  #
  # Assumes that range_name is also a column in bg_df (created by make_bg_booleans)
  
  if (! range_name %in% colnames(bg_df)) {
    abort(str_c(range_name, " column not contained in argument `bg_df` of ",
                            "`calculate_percent_in_given_range`."))
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
    c(str_c("percent_", range_name, "_full_day"),
      str_c("percent_", range_name, "_nighttime"),
      str_c("percent_", range_name, "_daytime"))
  
  percent_in_range_features
}

################################## STATS FOR GIVEN RANGES ######################

calculate_stat_in_given_range <- function(bg_df, range_name, f, stat_name) {
  # Applies f to valid BG readings that fall in the supplied
  # range_name (e.g. "very_high") for each day interval
  #
  # Assumes that range_name is also a column in bg_df (created by make_bg_booleans)
  
  if (! range_name %in% colnames(bg_df)) {
    abort(str_c(range_name, " column not contained in argument `bg_df` of ",
                "`calculate_percent_in_given_range`."))
  }
  
  if (! is_function(f)) {
    abort(str_c("Argument `f` of `calculate_stat_in_given_range` is not a ",
                "function"))
  }
  
  in_range <- pull(bg_df, range_name)
  nighttime <- pull(bg_df, "nighttime")
  bg <- pull(bg_df, "glucose")
  
  bg_in_range_full_day <- bg[in_range] 
  bg_in_range_nighttime <- bg[in_range & nighttime]
  bg_in_range_daytime <- bg[in_range & (! nighttime)]
  
  full_day_stat <- f(bg_in_range_full_day %>% check_if_sufficient_length)
  nighttime_stat <- f(bg_in_range_nighttime %>% check_if_sufficient_length)
  daytime_stat <- f(bg_in_range_daytime %>% check_if_sufficient_length)
  
  stat_in_range_features <- 
    c(full_day_stat,
      nighttime_stat,
      daytime_stat)
  
  names(stat_in_range_features) <- 
    c(str_c(range_name, "_", stat_name, "_bg_full_day"),
      str_c(range_name, "_", stat_name, "_bg_nighttime"),
      str_c(range_name, "_", stat_name, "_bg_daytime"))
  
  stat_in_range_features
}

################################## DEPRECATED ##################################

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


################################## AUC #########################################

# calculate_auc <- function(bg_df) {
#   # Calculates the AUC for a given day of CGM data, and scales the AUC by the
#   # number of valid recordings for that day
#   
#   time_diffs <- bg_df$datetime %>% 
#                   diff() %>% 
#                   time_length(unit = "minutes") %>% 
#                   abs() %>% 
#                   round()
#   
#   minutes_per_record <- calculate_minutes_per_record(bg_df)
#   
#   max_records_per_day <- 24 * calculate_records_per_hour(minutes_per_record)
#   
#   glucose <- bg_df$glucose
#   
#   # if time_diff is the same as normal recording interval, calculate auc, else 0
#   aucs <- pmap(list(head(glucose, -1),
#                     glucose[-1],
#                     time_diffs),
#        ~ if (..3 == minutes_per_record) {..3 * (..1 + ..2) / 2} else 0)
#   
#   aucs <- unlist(aucs)
#   
#   n_valid_aucs <- sum(aucs != 0)
#   
#   # if 288 possible recordings in a day, 287 possible aucs can be recorded
#   sum(aucs) * ((max_records_per_day - 1)  / n_valid_aucs)
# }
# 
# calculate_average_auc <- function(bg_df, range_name = NULL) {
#   
#   # if (! is.null(range_name)) {
#   #   
#   #   if (! range_name %in% colnames(bg_df)) {
#   #     abort(str_c(range_name, " column not contained in argument `bg_df` of `calculate_average_auc`."))
#   #   }
#   #   
#   #   in_range <- pull(bg_df, range_name)
#   #   
#   #   minutes_per_record <- calculate_minutes_per_record(bg_df)
#   #   
#   #   records_per_half_hour <- ceiling(60 / minutes_per_record)
#   #   
#   #   # collapse in_range 0 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 into 0 2 000 3 000 4
#   #   # if the numbers are >= records_per_half_hour, expand into corresponding number of 1s
#   #   # else expand into corresponding numbers of 0s
#   #
#   #   in_range_at_least_half_hour <- 
#   #   
#   #   bg_df <-
#   #     filter(bg_df, in_range_at_least_half_hour)
#   # }
#   
#   average_auc <- 
#     bg_df %>% 
#     group_by(cgm_day) %>% 
#     group_map(~ calculate_auc(.x)) %>% 
#     unlist() %>% 
#     mean()
#   
#   names(average_auc) <- str_c("auc_", range_name)
#   
#   average_auc
# }