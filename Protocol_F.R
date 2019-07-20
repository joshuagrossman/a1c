source("lib/source.R")

dfs <- make_dfs("data/Protocol_F/Data Tables")

#View(dfs)

cleaned_dfs <- list()

cleaned_dfs$cgm <-
  dfs$FDataCGM %>% 
  mutate(!is.na(Glucose)) %>% 
  select(RecID, PtID, DeviceDaysFromEnroll, DeviceTm, Glucose)

FBaseline_cleaned <- 
  dfs$FBaseline %>% 
  select(PtID, Gender, Weight, WeightUnits, Height, HeightUnits)

FPtRoster_cleaned <-
  dfs$FPtRoster %>% 
  select(PtID, RaceProtF)

FFinalVisit_cleaned <- 
  dfs$FFinalVisit %>% 
  mutate(a1c = as.double(HbA1c)) %>% 
  select(PtID, a1c, days_since_enrollment = TestDaysFromEnroll)

FSampleResults_cleaned <-
  dfs$FSampleResults %>% 
  filter(ResultName == "HbA1c") %>% 
  mutate(a1c = as.double(Value)) %>% 
  select(PtID, a1c, days_since_enrollment = CollectionDaysFromEnroll)

cleaned_dfs$measurements <-
  full_join(FBaseline_cleaned, FPtRoster_cleaned, by = "PtID")

a1c <- bind_rows(FFinalVisit_cleaned, FSampleResults_cleaned)

cleaned_dfs$a1c <- a1c

map2(cleaned_dfs,
     names(cleaned_dfs),
     ~ write_csv(.x, str_c("data/Protocol_F/clean/", .y, ".csv")))
