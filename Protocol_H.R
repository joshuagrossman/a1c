################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

DATA_PATH <- "data/Protocol_H/Data Tables"
CLEANED_DATA_PATH <- "data/Protocol_H/clean"
DATA_NAME <- "Protocol_H"
STUDY_START_DATE <- "03-01-2015"

CGM_NAME <- "Dexcom G4"

################################# READ IN DATA #################################

dfs <- make_dfs(DATA_PATH)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dfs$HDeviceCGM %>% 
  filter(!is.na(GlucoseValue)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDtTmDaysFromEnroll, 
                                        DeviceTm,
                                        STUDY_START_DATE)) %>% 
  select(id = PtID, 
         datetime, 
         glucose = GlucoseValue)

################################# HBA1C DATA ###################################

HLocalHbA1c_cleaned <- 
  dfs$HLocalHbA1c %>% 
  mutate(a1c = as.double(HbA1cTestRes),
         date = convert_to_date(HbA1cTestDtDaysAfterEnroll,
                                STUDY_START_DATE)) %>% 
  select(id = PtID, 
         a1c, 
         date)

Sample_cleaned <- 
  dfs$Sample %>% 
  filter(Analyte == "HBA1C") %>% 
  mutate(a1c = as.double(Value),
         date = convert_to_date(CollectionDtDaysFromEnroll,
                                STUDY_START_DATE)) %>% 
  select(id = PtID, 
         a1c, 
         date)

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
  select(id = PtID, 
         age = AgeAsOfEnrollDt)

HScreening_cleaned <-
  dfs$HScreening %>% 
  select(id = PtID, 
         gender = Gender, 
         ethnicity = Ethnicity, 
         race = Race, 
         weight = Weight, 
         height = Height)

cleaned_dfs$measurements <-
  full_join(HPtRoster_cleaned, HScreening_cleaned, by = "id") %>% 
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
