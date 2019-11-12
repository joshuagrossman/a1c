# Loads necessary libraries and data for CGM analysis

# NOTE: This takes awhile to run. Run once, and use the resulting output
# `data/features.csv` in any downstream analysis.

source("lib/feature_extraction/cgm_feature_extraction.R")
source("lib/feature_extraction/merging_cgm_data.R")

options(digits = 3, scipen = 999)

datasets <- c(
  "RT_CGM_Randomized_Clinical_Trial",
  "CMetforminDataset",
  "FLEX",
  "Protocol_F",
  "Protocol_H"
)

import_data <- function(dataset_name) {
  # Input: Dataset name
  
  # Output: List containing $cgm, $a1c, and $measurements data for the dataset.
  # See `load_file` definition for data format of each list item.
  
  cgm <- load_file(str_c("data/", dataset_name, "/clean/cgm.csv"), 
                   CGM_HEADERS)
  
  a1c <- load_file(str_c("data/", dataset_name, "/clean/a1c.csv"), 
                   A1C_HEADERS)
  
  measurements <- load_file(str_c("data/", dataset_name, 
                                  "/clean/measurements.csv"), 
                            NULL)
  
  list(cgm = cgm, a1c = a1c, measurements = measurements)
}

create_feature_df <- function(dataset_name) {
  # Input: Name of the dataset to import
  
  # Output: A dataframe with each row representing a single patient's a1c 
  # measurement with associated cgm and demographic data. Patients may have
  # more than one a1c measurement.
  
  data_list <- import_data(dataset_name)
  cgm <- data_list$cgm
  a1c <- data_list$a1c
  measurements <- data_list$measurements
  
  # creates list of merged cgm data with a1c and measurement data
  combined <- 
    combine_cgm_data_with_other_data(cgm, a1c, other_data_name = "a1c") %>% 
    split_patients_by_valid_a1c() %>% 
    combine_cgm_data_with_other_data(measurements, 
                                     other_data_name = "measurements")
  
  # creates feature_df for each list merged data
  combined_with_features <- 
    map(combined, 
        ~ append(.x, list(cgm_features = 
                            make_cgm_feature_df(.x, data_name = "cgm_data"))))
  
  # cleans and rbinds list of merged data
  feature_df_as_list <- 
    map(combined_with_features, 
        ~ c(id = str_c(.[["id"]], "_", dataset_name),
            a1c_date = .[["a1c_date"]],
            a1c_value = .[["a1c_value"]],
            .[["cgm_features"]], 
            .[["measurements"]]))
  
  print(paste0("Data from ", dataset, " has been imported."))
  reduce(feature_df_as_list, bind_rows)
}

bind_feature_df <- function(datasets) {
  # Input: Datasets to include in final feature_df
  
  # Output: Final feature_df, with each row corresponding to a single a1c
  # measurement for a single patient
  
  feature_df_as_list <- future_map(datasets, create_feature_df)

  feature_df <- reduce(feature_df_as_list, bind_rows)

  write_csv(feature_df, "data/features.csv")
}

bind_feature_df(datasets)