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
  
  if (! all(colnames(df) %in% required_headers)) {
    abort(str_c("Datafile at `filepath` must contain ", required_headers))
  }
  
  group_by(df, id) %>% 
    group_map(~ list(id = .y, data = .x))
}
