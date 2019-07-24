################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/BSevereHypoDataset/Data Tables"
cleaned_data_path <- "data/BSevereHypoDataset/clean"
data_name <- "BSevereHypoDataset"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <- 
  dfs$BDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDaysFromEnroll, DeviceTm)) %>% 
  select(RecID, PtID, datetime, Glucose)

################################# HBA1C DATA ###################################

cleaned_dfs$a1c <- 
  dfs$BSampleResults %>% 
  filter(ResultName == "HbA1c") %>% 
  mutate(days_since_enrollment = convert_to_date(0)) %>% 
  select(PtID, a1c = Value, days_since_enrollment)

################################ STATIC MEASUREMENTS ###########################

BDemoLifeDiabHxMgmt_cleaned <-
  dfs$BDemoLifeDiabHxMgmt %>% 
  select(PtID, Gender, Ethnicity, Race)

BMedChart_cleaned <- 
  dfs$BMedChart %>% 
  select(PtID, Weight, WeightUnits, Height, HeightUnits)

BPtRoster_cleaned <- 
  dfs$BPtRoster %>% 
  select(PtID, BCaseControlStatus)

cleaned_dfs$measurements <-
  full_join(BDemoLifeDiabHxMgmt_cleaned, BMedChart_cleaned, by = "PtID") %>% 
  full_join(BPtRoster_cleaned, by = "PtID")

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)
     
