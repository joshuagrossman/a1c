################################## LIBRARIES ###################################

source("lib/cleaning_data/clean.R")

################################## GLOBALS #####################################

DATA_PATH <- "data/CMetforminDataset/Data Tables"
CLEANED_DATA_PATH <- "data/CMetforminDataset/clean"
DATA_NAME <- "CMetforminDataset"
STUDY_START_DATE <- "07-02-2014"

################################# READ IN DATA #################################

dfs <- make_dfs(DATA_PATH)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <- 
  dfs$CDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDaysFromEnroll, 
                                        DeviceTm, 
                                        STUDY_START_DATE)) %>% 
  select(id = PtID, 
         datetime, 
         glucose = Glucose)

################################ HEIGHT / WEIGHT ##############################

# CRandomization_hw <-
#   dfs$CRandomization %>% 
#   select(id = PtID, 
#          visit = NA,
#          weight = Weight, 
#          height = Height)

CFollowUpVisit_hw <-
  dfs$CFollowUpVisit %>% 
  select(id = PtID,
         visit = Visit,
         weight = Weight,
         height = Height)

# CPostTreatmentVisit_hw <- 
#   dfs$CPostTreatmentVisit %>% 
#   select(id = PtID,
#          visit = NA,
#          height = Height,
#          weight = Weight)

CScreening_hw <- 
  dfs$CScreening %>% 
  # no Visit column, so manually entering this string
  mutate(visit = "Screening") %>% 
  select(id = PtID,
         visit,
         height = Height,
         weight = Weight)

hw_cleaned <-
  bind_rows(CFollowUpVisit_hw, CScreening_hw)

###################### HBA1C DATA W/ HEIGHT/WEIGHT ###########################

CLabHbA1c_cleaned <- 
  dfs$CLabHbA1c %>% 
  mutate(date = convert_to_date(TestDaysFromEnroll,
                                STUDY_START_DATE)) %>% 
  select(id = PtID,
         a1c = HbA1c,
         date,
         visit = Visit)

# Impossible to determine when these a1c measurements occurred, seems like they
# are before the study start date
# SampleResults_cleaned <- 
#   dfs$SampleResults %>% 
#   filter(ResultName == "HbA1c") %>% 
#   select(id, a1c = Value)

cleaned_dfs$a1c <- 
  CLabHbA1c_cleaned %>% 
  left_join(hw_cleaned, by = c("id", "visit"))

############################ OTHER MEASUREMENTS ################################

CScreening_cleaned <-
  dfs$CScreening %>% 
  mutate(diag_date = str_replace(T1DDiagYrMon, 
                                 "(\\d{4})(\\d{2})", 
                                 "\\2-01-\\1") %>% 
                     mdy(),
         estimated_dob = diag_date - years(DiagT1DAge),
         age = interval(estimated_dob, 
                                  mdy(STUDY_START_DATE)) / years(1),
         years_since_diag = age - as.double(DiagT1DAge)) %>% 
  select(id = PtID, 
         gender = Gender, 
         ethnicity = Ethnicity, 
         race = Race, 
         age,
         years_since_diag)

CTrtGroupUnmasked_cleaned <-
  dfs$CTrtGroupUnmasked %>% 
  select(id = PtID, 
         metformin = TrtGroup)

cgm_name <-
  dfs$CDataCGM %>% 
  group_by(PtID) %>% 
  select(id = PtID, 
         cgm_name = CCGMDeviceType) %>% 
  summarize(cgm_name = unique(cgm_name))

cleaned_dfs$measurements <- 
  CScreening_cleaned %>%  
  full_join(CTrtGroupUnmasked_cleaned, by = "id") %>% 
  full_join(cgm_name, by = "id") %>% 
  mutate(datafile = str_c(.$id, ".csv"),
         study_start_date = STUDY_START_DATE,
         dataset = DATA_NAME)



################################ DATA CHECKS ###################################

only_unique_patients <- nrow(cleaned_dfs$measurements) == 
  length(unique(cleaned_dfs$measurements$id))
  
if (!only_unique_patients) {
  warning("Non-unique patient data.")
}

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, CLEANED_DATA_PATH)

# split_cgm_into_patients(CLEANED_DATA_PATH)

