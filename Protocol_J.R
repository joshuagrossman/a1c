################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/Protocol_J/Data Tables"
cleaned_data_path <- "data/Protocol_Jclean"
data_name <- "Protocol_J"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dfs$IDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDtDaysFromEnroll, DeviceTm)) %>% 
  select(PtID, datetime, Glucose)

################################# HBA1C DATA ###################################

cleaned_dfs$a1c <-
  dfs$IScreening %>% 
  mutate(date = convert_to_date(TestDtDaysAfterEnroll)) %>% 
  select(PtID, date, HbA1c)

################################# INSULIN DATA #################################

# TODO(Josh): Fix this if insulin data is useful

# cleaned_dfs$insulin <-
#   dfs$HDataPump %>%
#   mutate(datetime = convert_to_datetime(DeviceDtDaysFromEnroll, DeviceTm)) %>%
#   select(PtID, datetime, ImmediateVol, ExtendedVol)

################################ STATIC MEASUREMENTS ###########################

IPtRoster_cleaned <-
  dfs$IPtRoster %>% 
  select(PtID, AgeAsOfEnrollDt)

IScreening_cleaned <-
  dfs$IScreening %>% 
  select(PtID, Gender, Ethnicity, Race, Weight, Height)

cleaned_dfs$measurements <-
  full_join(IPtRoster_cleaned, IScreening_cleaned, by = "PtID")

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)
