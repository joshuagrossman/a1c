################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/Protocol_F/Data Tables"
cleaned_data_path <- "data/Protocol_F/clean"
data_name <- "Protocol_F"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dfs$FDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDaysFromEnroll, DeviceTm)) %>% 
  select(RecID, PtID, datetime, Glucose)

################################# HBA1C DATA ###################################

FFinalVisit_cleaned <- 
  dfs$FFinalVisit %>% 
  mutate(a1c = as.double(HbA1c),
         date = convert_to_date(TestDaysFromEnroll)) %>% 
  select(PtID, a1c, date)

FSampleResults_cleaned <-
  dfs$FSampleResults %>% 
  filter(ResultName == "HbA1c") %>% 
  mutate(a1c = as.double(Value),
         date = convert_to_date(CollectionDaysFromEnroll)) %>% 
  select(PtID, a1c, date)

cleaned_dfs$a1c <- bind_rows(FFinalVisit_cleaned, FSampleResults_cleaned)

################################ STATIC MEASUREMENTS ###########################

FBaseline_cleaned <- 
  dfs$FBaseline %>% 
  select(PtID, Gender, Weight, WeightUnits, Height, HeightUnits)

FPtRoster_cleaned <-
  dfs$FPtRoster %>% 
  select(PtID, RaceProtF)

cleaned_dfs$measurements <-
  full_join(FBaseline_cleaned, FPtRoster_cleaned, by = "PtID")

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)
