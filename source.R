################################## LIBRARIES ###################################

library(tidyverse)
library(fs)
library(lubridate)

################################## READING FILES ###############################

is_csv <- function(file) {
  first_line <- read_lines(file, n_max = 1)
  n_chars <- str_length(first_line)
  n_commas <- str_count(first_line, ",")
  n_pipes <- str_count(first_line, "\\|")
  
  return((n_commas/n_chars) > (n_pipes/n_chars))
}

parse_file <- function(filename) {
  if (is_csv(filename)) {
    return(read_csv(filename, cols(.default = "c")))
  }
  
  read_delim(filename, "|")
}

make_dfs <- function(path) {
  filenames <- dir_ls(path)
  dfs <- map(filenames, parse_file)
  base_filenames <- str_replace(filenames, ".*/(.+)\\.txt", "\\1")
  names(dfs) <- base_filenames
  dfs
}

################################## MUTATING CSVS ###############################

convert_to_date <- function(days_since_enrollment) {
  mdy("1/1/2019") + days(days_since_enrollment)
}

convert_to_datetime <- function(days_since_enrollment, device_time) {
  mdy_hms(str_c("1/1/2019 ", device_time)) + days(days_since_enrollment)
}

################################## WRITING CSVS ################################

write_cgm_csv_for_patient <- function(df, dataset_name) {
  patient_id <- df$patient_id[1]
  write_csv(df, str_c("data/", 
                      dataset_name, 
                      "/clean/cgm_by_patient/", 
                      patient_id, 
                      ".csv"))
}

split_into_csvs <- function(cgm_data_path, dataset_name) {
  cgm <- 
    read_csv(cgm_data_path) %>% 
    mutate(dataset = dataset_name,
           patient_id = PtID) %>% 
    group_by(PtID)
  
  group_walk(cgm, ~ write_cgm_csv_for_patient(., dataset_name))
}