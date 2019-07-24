################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/Protocol_H/Data Tables"
cleaned_data_path <- "data/Protocol_H/clean"
data_name <- "Protocol_H"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dfs$HDeviceCGM %>% 
  filter(!is.na(GlucoseValue)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDtTmDaysFromEnroll, DeviceTm)) %>% 
  select(RecID, PtID, datetime, GlucoseValue)

################################# HBA1C DATA ###################################

HLocalHbA1c_cleaned <- 
  dfs$HLocalHbA1c %>% 
  mutate(a1c = as.double(HbA1cTestRes),
         date = convert_to_date(HbA1cTestDtDaysAfterEnroll)) %>% 
  select(PtID, a1c, date)

Sample_cleaned <- 
  dfs$Sample %>% 
  filter(Analyte == "HBA1C") %>% 
  mutate(a1c = as.double(Value),
         date = convert_to_date(CollectionDtDaysFromEnroll)) %>% 
  select(PtID, a1c, date)

cleaned_dfs$a1c <- bind_rows(HLocalHbA1c_cleaned, Sample_cleaned)

################################# INSULIN DATA #################################

# TODO(Josh): Fix this if insulin data is useful

# cleaned_dfs$insulin <-
#   dfs$HDeviceBolus %>% 
#   mutate(datetime = convert_to_datetime(DeviceDtTmDaysFromEnroll, DeviceTm)) %>% 
#   select(PtID, datetime, InjValue)

################################ STATIC MEASUREMENTS ###########################

HPtRoster_cleaned <- 
  dfs$HPtRoster %>% 
  select(PtID, AgeAsOfEnrollDt)

HScreening_cleaned <-
  dfs$HScreening %>% 
  select(PtID, Gender, Ethnicity, Race, Weight, Height)

cleaned_dfs$measurements <-
  full_join(HPtRoster_cleaned, HScreening_cleaned, by = "PtID")

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)
