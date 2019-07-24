################################## LIBRARIES ###################################

source("lib/source.R")

################################## GLOBALS #####################################

data_path <- "data/CMetforminDataset/Data Tables"
cleaned_data_path <- "data/CMetforminDataset/clean"
data_name <- "CMetforminDataset"

################################# READ IN DATA #################################

dfs <- make_dfs(data_path)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <- 
  dfs$CDataCGM %>% 
  filter(!is.na(Glucose)) %>% 
  mutate(datetime = convert_to_datetime(DeviceDaysFromEnroll, DeviceTm)) %>% 
  select(PtID, datetime, Glucose)

################################# HBA1C DATA ###################################

cleaned_dfs$a1c <- 
  dfs$CLabHbA1c %>% 
  mutate(date = convert_to_date(TestDaysFromEnroll)) %>% 
  select(PtID, date, HbA1c)

# Impossible to determine when these a1c measurements occurred, seems like they
# are before the study start date
# SampleResults_cleaned <- 
#   dfs$SampleResults %>% 
#   filter(ResultName == "HbA1c") %>% 
#   select(PtID, a1c = Value)

################################ STATIC MEASUREMENTS ###########################
  
CRandomization_cleaned <-
  dfs$CRandomization %>% 
  select(PtID, Weight, Height)

CScreening_cleaned <-
  dfs$CScreening %>% 
  mutate(diag_date = str_replace(T1DDiagYrMon, 
                                 "(\\d{4})(\\d{2})", 
                                 "\\2-01-\\1") %>% 
                     mdy(),
         estimated_dob = diag_date - years(DiagT1DAge),
         # study occurred end of 2013 / beginning of 2014
         estimated_age = interval(estimated_dob, 
                                  mdy("1/1/2014")) / years(1)) %>% 
  select(PtID, Gender, Ethnicity, Race, estimated_age)

CTrtGroupUnmasked_cleaned <-
  dfs$CTrtGroupUnmasked %>% 
  select(PtID, TrtGroup)

cleaned_dfs$measurements <- 
  left_join(CRandomization_cleaned, CScreening_cleaned, by = "PtID") %>% 
  full_join(CTrtGroupUnmasked_cleaned, by = "PtID") %>% 
  mutate(datafile = str_c(.$PtID, ".csv"))

only_unique_patients <- nrow(cleaned_dfs$measurements) == length(unique(cleaned_dfs$measurements$PtID))
  
if (!only_unique_patients) {
  warning("Non-unique patient data.")
}

################################ WRITE CSVS TO FILE ############################

write_main_csvs(cleaned_dfs, cleaned_data_path)

split_into_csvs(str_c(cleaned_data_path, "/cgm.csv"), data_name)

