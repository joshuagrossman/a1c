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

################################## GLOBALS #####################################

MAX_RECORDS_PER_DAY_5_MINUTE <- 288
  
MAX_RECORDS_PER_DAY_15_MINUTE <- 96

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

