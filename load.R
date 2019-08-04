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
  
  if(is_dir(filename)) {
    warning(str_c("Can't parse a directory: ", filename))
    return(NULL)
  }
  
  if (is_csv(filename)) {
    return(read_csv(filename, col_types = cols(.default = "c")))
  }
  
  read_delim(filename, "|", col_types = cols(.default = "c"))
}

make_dfs <- function(path) {
  filenames <- dir_ls(path)
  dfs <- map(filenames, parse_file)
  base_filenames <- str_replace(filenames, ".*/(.+)\\.[:alpha:]{3}", "\\1")
  names(dfs) <- base_filenames
  dfs
}

################################## MUTATING CSVS ###############################

convert_to_date <- function(days_since_enrollment, study_start_date) {
  mdy(study_start_date) + days(days_since_enrollment)
}

convert_to_datetime <- function(days_since_enrollment, device_time, study_start_date) {
  mdy_hms(str_c(study_start_date, " ", device_time)) + days(days_since_enrollment)
}

inches_to_cm <- function(inches) {
  2.54 * as.double(inches)
}

lbs_to_kg <- function(lbs) {
  0.4536 * as.double(lbs)
}

################################## WRITING CSVS ################################

write_main_csvs <- function(cleaned_dfs, cleaned_data_path) {
  map2(cleaned_dfs,
       names(cleaned_dfs),
       ~ write_csv(.x, str_c(cleaned_data_path, "/", .y, ".csv")))
}

write_cgm_csv_for_patient <- function(df, cleaned_data_path) {
  patient_id <- df$id[1]
  df_no_id <- select(df, -id)
  write_csv(df_no_id, str_c(cleaned_data_path,
                      "/cgm_by_patient/",
                      patient_id, 
                      ".csv"))
}

split_cgm_into_patients <- function(cleaned_data_path) {
  cgm <- 
    read_csv(str_c(cleaned_data_path, "/cgm.csv")) %>% 
    group_by(id)
  
  group_walk(cgm, 
             ~ write_cgm_csv_for_patient(., cleaned_data_path), 
             keep = T)
}