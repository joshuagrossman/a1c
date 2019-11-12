################################## LIBRARIES ###################################

source("lib/cleaning_data/clean.R")

################################## GLOBALS #####################################

DATA_PATH <- "data/RT_CGM_Randomized_Clinical_Trial/DataTables"
CLEANED_DATA_PATH <- "data/RT_CGM_Randomized_Clinical_Trial/clean"
DATA_NAME <- "RT_CGM_Randomized_Clinical_Trial"
STUDY_START_DATE <- "12-01-2006"

CGM_NAME <- "Freestyle Navigator or DexCom STS or Medtronic Paradigm/Guardian"

################################# READ IN DATA #################################

dfs <- make_dfs(DATA_PATH)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <-
  dir_map(str_c(DATA_PATH, "/cgm"), read_csv) %>% 
  reduce(bind_rows) %>% 
  select(id = PtID,
         datetime = DeviceDtTm,
         glucose = Glucose)

################################# HEIGHT / WEIGHT ##############################

tblAPtSummary_hw <-
  dfs$tblAPtSummary %>% 
  transmute(id = PtID,
         # not clear when this height/weight was recorded, potentially at randomization
         visit = NA_character_,
         height = Height, 
         weight = Weight)

# this df is missing from the data, but included in the study description
# HScreening_hw <-
#   dfs$HScreening %>% 
#   transmute(id = PtID,
#          visit = "Screening",
#          height = Height,
#          weight = Weight)

# cleaned_hw <- bind_rows(tblAPtSummary_hw, HScreening_hw)


################################# HBA1C DATA ###################################

cleaned_dfs$a1c <-
  dfs$tblALabHbA1c %>% 
  select(id = PtID, 
         a1c = LabA1cResult,
         date = LabHbA1cDt,
         visit = Visit) %>%
  left_join(tblAPtSummary_hw, by = c("id", "visit"))

################################# INSULIN DATA #################################

# TODO(Josh): Fix this if insulin data is useful

# cleaned_dfs$insulin <-
#   dfs$HDeviceBolus %>% 
#   mutate(datetime = convert_to_datetime(DeviceDtTmDaysFromEnroll, DeviceTm)) %>% 
#   select(PtID, datetime, InjValue)

################################ STATIC MEASUREMENTS ###########################

cleaned_dfs$measurements <-
  dfs$tblAPtSummary %>% 
  select(id = PtID, 
         gender = Gender, 
         age = AgeAsOfRandDt, 
         race = Race, 
         ethnicity = Ethnicity) %>% 
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

# split_cgm_into_patients(CLEANED_DATA_PATH)
