################################## LIBRARIES ###################################

library(tidyverse)
library(fs)
library(rlang)
library(lubridate)
library(furrr)

plan(multisession(workers = availableCores()))

################################## GLOBALS #####################################

CGM_HEADERS <- c("id", "glucose", "datetime")

A1C_HEADERS <- c("id", "a1c", "date")

################################## LOAD DATA ###################################

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
    if (! all(colnames(df) %in% required_headers)) {
      abort(str_c("Datafile at `filepath` must contain ", required_headers))
    }
  }
  
  group_by(df, id) %>% 
    group_map(~ list(id = pull(.y, "id"), data = .x))
}



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