source("lib/source.R")

dfs <- make_dfs("data/BSevereHypoDataset/Data Tables")

#View(dfs)

cleaned_dfs <- list()

cleaned_dfs$BDataCGM <- 
  dfs$BDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  # group_by(PtID) %>% 
  # mutate(Day = DeviceDaysFromEnroll - min(DeviceDaysFromEnroll)) %>% 
  select(RecID, PtID, DeviceDaysFromEnroll, DeviceTm, Glucose)

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

cleaned_dfs$a1c <- 
  dfs$BSampleResults %>% 
  filter(ResultName == "HbA1c") %>% 
  mutate(days_since_enrollment = 0) %>% 
  select(PtID, a1c = Value, days_since_enrollment)

map2(cleaned_dfs,
     names(cleaned_dfs),
     ~ write_csv(.x, str_c("data/BSevereHypoDataset/clean/", .y, ".csv")))
     
