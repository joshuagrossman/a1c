################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

DATA_PATH <- "data/Protocol_F/Data Tables"
CLEANED_DATA_PATH <- "data/Protocol_F/clean"
DATA_NAME <- "Protocol_F"
STUDY_START_DATE <- "01-01-2015"

CGM_NAME <- "Freestyle Libre Pro Flash"

################################# READ IN DATA #################################

dfs <- make_dfs(DATA_PATH)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dfs$FDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDaysFromEnroll, 
                                        DeviceTm,
                                        STUDY_START_DATE)) %>% 
  select(id = PtID, 
         datetime, 
         glucose = Glucose)

################################# HBA1C DATA ###################################

FFinalVisit_cleaned <- 
  dfs$FFinalVisit %>% 
  mutate(a1c = as.double(HbA1c),
         date = convert_to_date(TestDaysFromEnroll, 
                                STUDY_START_DATE)) %>% 
  select(id = PtID, 
         a1c, 
         date)

FSampleResults_cleaned <-
  dfs$FSampleResults %>% 
  filter(ResultName == "HbA1c") %>% 
  mutate(a1c = as.double(Value),
         date = convert_to_date(CollectionDaysFromEnroll,
                                STUDY_START_DATE)) %>% 
  select(id = PtID, 
         a1c, 
         date)

cleaned_dfs$a1c <- bind_rows(FFinalVisit_cleaned, FSampleResults_cleaned)

################################ STATIC MEASUREMENTS ###########################

FBaseline_cleaned <- 
  dfs$FBaseline %>% 
  transmute(id = PtID, 
            gender = Gender, 
            weight = ifelse(WeightUnits == "lbs", lbs_to_kg(Weight), Weight),
            height = ifelse(HeightUnits == "in", inches_to_cm(Height), Height))

FPtRoster_cleaned <-
  dfs$FPtRoster %>% 
  mutate(ethnicity = str_replace(RaceProtF, "(.+) (.+)", "\\1"),
         race = str_replace(RaceProtF, "(.+) (.+)", "\\2")) %>% 
  select(id = PtID, 
         race,
         ethnicity)

cleaned_dfs$measurements <-
  full_join(FBaseline_cleaned, FPtRoster_cleaned, by = "id") %>% 
  mutate(datafile = str_c(.$id, ".csv"),
         study_start_date = STUDY_START_DATE,
         dataset = DATA_NAME,
         cgm_name = CGM_NAME)

################################ DATA CHECKS ###################################

only_unique_patients <- nrow(cleaned_dfs$measurements) == 
  length(unique(cleaned_dfs$measurements$id))

if (!only_unique_patients) {
  warning("Non-unique patient data.")
}

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, CLEANED_DATA_PATH)

split_cgm_into_patients(CLEANED_DATA_PATH)
