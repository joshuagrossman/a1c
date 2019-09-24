################################## LIBRARIES ###################################

library(tidyverse)
library(fs)
library(rlang)
library(lubridate)
library(furrr)

# Use parallelization responsibly
plan(multicore(workers = availableCores() %/% 2 + 1))

################################## GLOBALS #####################################

CGM_HEADERS <- c("id", "glucose", "datetime")

A1C_HEADERS <- c("id", "a1c", "date")

################################## LOAD DATA ###################################

make_id_list <- function(df) {
  # Input: Data with an `id` column
  #
  # Output: A list of sublists, each containing an id and its associated data
  
  group_by(df, id) %>% 
    group_map(~ list(id = pull(.y, "id"), data = .x))
}

load_file <- function(filepath, required_headers) {
  # Input: Filepath to a CSV containing data for more than one patient. The
  # CSV must contain required_headers as columns.
  #
  # Output: A list of lists, with each list containing the patient_id and all
  # associated data as a df.
  
  if (path_ext(filepath) != "csv") {
    abort("Datafile at `filepath` must be a CSV.")
  }
  
  df <- read_csv(filepath)
  
  if (! is.null(required_headers)) {
    if (! all(required_headers %in% colnames(df))) {
      abort(str_c("Datafile at `filepath` must contain ", required_headers))
    }
  }
  
  make_id_list(df)
}

load_single_patient_file <- function(filepath) {

  if (path_ext(filepath) != "csv") {
    abort("CGM datafile must be a CSV.")
  }

  df <- read_csv(filepath)
  
  if (! is.null(required_headers)) {
    if (! all(colnames(df) %in% required_headers)) {
      abort(str_c("Datafile at `filepath` must contain ", required_headers))
    }
  }

  id <- path_ext_remove(path_file(filepath))

  list(id = id, data = df)
}

load_dir_of_patients <- function(dirpath) {

  if (! is_dir(dirpath)) {
    abort("CGM data folder must be a directory.")
  }

  dir_map(dirpath, load_single_patient_file)
}