################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/RT_CGM_Randomized_Clinical_Trial/DataTables"
cleaned_data_path <- "data/RT_CGM_Randomized_Clinical_Trial/clean"
data_name <- "RT_CGM_Randomized_Clinical_Trial"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dir_map(str_c(data_path, "/cgm"), read_csv) %>% 
  reduce(bind_rows)

################################# HBA1C DATA ###################################

cleaned_dfs$a1c <-
  dfs$tblALabHbA1c %>% 
  select(PtID, LabHbA1cDt, LabA1cResult)

################################# INSULIN DATA #################################

# TODO(Josh): Fix this if insulin data is useful

# cleaned_dfs$insulin <-
#   dfs$HDeviceBolus %>% 
#   mutate(datetime = convert_to_datetime(DeviceDtTmDaysFromEnroll, DeviceTm)) %>% 
#   select(PtID, datetime, InjValue)

################################ STATIC MEASUREMENTS ###########################

cleaned_dfs$measurements <-
  dfs$tblAPtSummary %>% 
  select(PtID, Gender, AgeAsOfRandDt, Race, Ethnicity, Height, Weight)

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)
