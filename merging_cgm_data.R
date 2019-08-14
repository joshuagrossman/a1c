source("lib/load.R")

################# COMBINING OTHER PATIENT DATA WITH CGM DATA ###################

create_list_with_name <- function(x, x_name) {
  l <- list(x)
  names(l) <- x_name
  l
}

combine_cgm_data_with_other_data <- function(cgm_data, 
                                     other_data,
                                     other_data_name,
                                     id_name = "id",
                                     data_name = "data") {
  # Input: Two lists of sub-lists, with each sublist containing the patient id
  # and associated cgm or other data.
  #
  # Output: One list of sub-lists, each containing id, other data, cgm data for a single patient.
  
  cgm_ids <- map_chr(cgm_data, ~ as.character(.[[id_name]]))
  other_ids <- map_chr(other_data, ~ as.character(.[[id_name]]))
  
  matching_id_indices <- match(cgm_ids, other_ids)
  
  filtered_cgm_data <- cgm_data
  filtered_cgm_data[is.na(matching_id_indices)] <- NULL
  
  matching_id_indices <- matching_id_indices[! is.na(matching_id_indices)]
  
  combined_data <- map2(.x = matching_id_indices, 
                        .y = filtered_cgm_data,
                        .f = ~ append(.y, 
                                      create_list_with_name(other_data[[.x]][[data_name]],
                                                            other_data_name)
                                      )
                        )
}

################################ A1C-SPECIFIC FILTERING ########################

identify_most_recent_a1c_with_cgm_data <- function(individual_patient_data,
                                                   data_name = "data",
                                                   a1c_name = "a1c",
                                                   date_name = "date",
                                                   datetime_name = "datetime") {
  
  # Identifies the date of the most recent a1c that falls inside the date range 
  # of cgm data
  #
  # TODO: Have this select the most recent a1c date with sufficient cgm data,
  # not just any amount of CGM data
  
  bg_df <- individual_patient_data[[data_name]]
  a1c <- individual_patient_data[[a1c_name]]
  
  max_cgm_date <- max(pull(bg_df, datetime_name))
  min_cgm_date <- min(pull(bg_df, datetime_name))
  
  most_recent_a1c <- a1c %>% 
    filter(between(.data[[date_name]], min_cgm_date, max_cgm_date)) %>% 
    arrange(desc(date)) %>% 
    slice(1)
  
  if (nrow(most_recent_a1c) == 0) {
    return(c(individual_patient_data,
                most_recent_a1c_value = NA, 
                most_recent_a1c_date = NA))
  }
  
  c(individual_patient_data, 
    most_recent_a1c_value = pull(most_recent_a1c, a1c_name),
    most_recent_a1c_date = pull(most_recent_a1c, date_name))
  
}

filter_cgm_by_date <- function(bg_df, end_date, max_days) {
  # Removes CGM data that falls after `end_date` or before `end_date` - `max_days`
    
  end_date_as_date <- as.POSIXct(end_date, origin="1970-01-01")
    
  bg_df %>%
    filter(between(datetime, 
                   end_date_as_date - days(max_days),
                   end_date_as_date - days(1)))
}

filter_individual_cgm_by_a1c <- function(individual_patient_data, max_days) {
  # Removes CGM data that falls after the most recent a1c date or before
  # `max_days` before the most recent a1c date
  
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