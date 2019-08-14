source("lib/load.R")

combine_a1c_and_cgm_data <- function(cgm_data, a1c_data) {
  # Input: Two lists of sub-lists, with each sublist containing the patient id
  # and associated cgm or a1c data.
  #
  # Output: One list of sub-lists, each containing id, a1c, cgm data for a single patient.
  
  cgm_ids <- map_chr(cgm_data, ~ as.character(.$id))
  a1c_ids <- map_chr(a1c_data, ~ as.character(.$id))
  
  matching_id_indices <- match(cgm_ids, a1c_ids)
  
  filtered_cgm_data <- cgm_data
  filtered_cgm_data[is.na(matching_id_indices)] <- NULL
  
  matching_id_indices <- matching_id_indices[! is.na(matching_id_indices)]
  
  combined_data <- map2(.x = matching_id_indices, 
                        .y = filtered_cgm_data,
                        .f = ~ c(.y, list(a1c = a1c_data[[.x]]$data)))
}

combine_measurements_and_cgm_data <- function(cgm_data, measurements) {
  # Input: Two lists of sub-lists, with each sublist containing the patient id
  # and associated cgm or a1c data.
  #
  # Output: One list of sub-lists, each containing id, a1c, cgm data for a single patient.
  
  cgm_ids <- map_chr(cgm_data, ~ as.character(.$id))
  measurements_ids <- map_chr(measurements, ~ as.character(.$id))
  
  matching_id_indices <- match(cgm_ids, measurements_ids)
  
  filtered_cgm_data <- cgm_data
  filtered_cgm_data[is.na(matching_id_indices)] <- NULL
  
  matching_id_indices <- matching_id_indices[! is.na(matching_id_indices)]
  
  combined_data <- map2(.x = matching_id_indices, 
                        .y = filtered_cgm_data,
                        .f = ~ c(.y, list(measurements = measurements[[.x]]$data)))
}

identify_most_recent_a1c_with_cgm_data <- function(individual_patient_data) {
  
  bg_df <- individual_patient_data$data
  a1c <- individual_patient_data$a1c
  
  max_cgm_date <- max(bg_df$datetime)
  min_cgm_date <- min(bg_df$datetime)
  
  most_recent_a1c <- a1c %>% 
    filter(between(date, min_cgm_date, max_cgm_date)) %>% 
    arrange(desc(date)) %>% 
    slice(1)
  
  if (nrow(most_recent_a1c) == 0) {
    return(list(individual_patient_data,
                most_recent_a1c_value = NA, 
                most_recent_a1c_date = NA))
  }
  
  c(individual_patient_data, 
    most_recent_a1c_value = most_recent_a1c$a1c,
    most_recent_a1c_date = most_recent_a1c$date)
  
}

filter_cgm_by_date <- function(bg_df, end_date, max_days) {
    
  end_date_as_date <- as.POSIXct(end_date, origin="1970-01-01")
    
  bg_df %>%
    filter(between(datetime, 
                   end_date_as_date - days(max_days),
                   end_date_as_date - days(1)))
}

filter_individual_cgm_by_a1c <- function(individual_patient_data, max_days) {
  
  recent_a1c_identifed <- identify_most_recent_a1c_with_cgm_data(individual_patient_data)
  
  most_recent_a1c_date <- recent_a1c_identifed$most_recent_a1c_date
    
  if (is.na(most_recent_a1c_date)) {
    return(c(recent_a1c_identifed, filtered_cgm_data = NA))
  }
  
  recent_a1c_identifed$filtered_cgm_data <- 
    filter_cgm_by_date(recent_a1c_identifed$data, 
                       most_recent_a1c_date,
                       max_days)
  
  recent_a1c_identifed
}

filter_all_cgm_by_a1c <- function(patient_data, max_days = 30) {
  
  map(patient_data, 
      filter_individual_cgm_by_a1c, 
      max_days)
}