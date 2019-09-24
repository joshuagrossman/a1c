source("lib/feature_extraction/load.R")

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
  # Input: Two lists of sub-lists created by load_file(), with each sublist 
  # containing the patient id and associated data.
  #
  # Output: One list of sub-lists, each containing patient id and all associated data.
  
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
                                                            other_data_name)))
  
  combined_data
}

################################ MERGING A1C DATA #############################

split_patients_by_valid_a1c <- function(individual_patient_data,
                                        id_name = "id",
                                        a1c_data_name = "a1c",
                                        cgm_data_name = "data",
                                        max_days = 30,
                                        extra_days = 7) {
  
  combined <- 
    map(individual_patient_data, 
        ~ extract_valid_a1cs_and_associated_cgm_data(.[[a1c_data_name]],
                                                   .[[cgm_data_name]],
                                                   .[[id_name]],
                                                   max_days,
                                                   extra_days)) %>% 
    flatten()
}

filter_cgm_by_date <- function(bg_df, end_date, max_days) {
  # Removes CGM data that falls after `end_date` or before `end_date` - `max_days`
    
  end_date_as_date <- as_date(end_date)
    
  bg_df %>%
    filter(between(as_date(datetime), 
                   end_date_as_date - days(max_days),
                   end_date_as_date - days(1)))
}

extract_valid_a1cs_and_associated_cgm_data <- function(a1c_df, 
                                                       bg_df,
                                                       id,
                                                       max_days,
                                                       extra_days = 7) {
  
  if (all(is.na(a1c_df))) {
    warning("No a1c data provided.")
    return(list(list(id = id,
                     a1c_value = NA,
                     a1c_date = NA,
                     cgm_data = bg_df)))
  }
  
  pmap(a1c_df, ~ list(id = id,
                      a1c_value = ..1,
                      a1c_date = ..2,
                      cgm_data = filter_cgm_by_date(bg_df,
                                                    ..2 + days(extra_days),
                                                    max_days)))
}


################################ DEPRECATED ####################################

# identify_a1cs_with_cgm_data <- function(a1c_df,
#                                         bg_df,
#                                         a1c_name = "a1c",
#                                         date_name = "date",
#                                         datetime_name = "datetime") {
#   
#   # Identifies the lab-measured a1cs that occurred at the same time BG was 
#   # measured using CGM
#   
#   cgm_datetimes <- pull(bg_df, datetime_name)
#     
#   max_cgm_date <- as_date(max(cgm_datetimes, na.rm = T))
#   min_cgm_date <- as_date(min(cgm_datetimes, na.rm = T))
#   
#   valid_a1cs <- a1c_df %>%
#     filter(! is.na(a1c)) %>% 
#     filter(between(as_date(date), min_cgm_date, max_cgm_date)) %>% 
#     arrange(desc(date))
#   
#   if (nrow(valid_a1cs) == 0) {
#     return(NA)
#   }
#   
#   tibble(
#     a1c_value = pull(valid_a1cs, a1c_name),
#     a1c_date = as_date(pull(valid_a1cs, date_name))
#   )
#   
# }