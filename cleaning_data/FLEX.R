################################## LIBRARIES ###################################

source("lib/cleaning_data/clean.R")

################################## GLOBALS #####################################

DATA_PATH <- "data/FLEX/Data Tables"
CLEANED_DATA_PATH <- "data/FLEX/clean"
DATA_NAME <- "FLEX"
STUDY_START_DATE <- "05-01-2014"

CGM_NAME <- "iPro2 Professional CGM"

################################# READ IN DATA #################################

dfs <- make_dfs(DATA_PATH)

#View(dfs)

cleaned_dfs <- list()

################################# CGM DATA #####################################

cleaned_dfs$cgm <- 
  dfs$raw_cgm %>% 
  mutate(glucose = as.numeric(Glucose),
         datetime = as.numeric(DeviceDtTm)) %>% 
  filter(!is.na(glucose)) %>% 
  mutate(datetime = as_datetime(datetime, origin = "1960-01-01")) %>% 
  select(id = newID, 
         datetime, 
         glucose)

################################ HEIGHT / WEIGHT ##############################

cleaned_hw <- 
  dfs$raw_measurements %>% 
  transmute(id = newID,
         bmi_z = as.numeric(bmiz),
         age = age_yrs,
         date = visit_date)

###################### HBA1C DATA W/ HEIGHT/WEIGHT ###########################

cleaned_dfs$a1c <-
  dfs$raw_a1c %>% 
  select(id = newID,
         date = visit_date,
         a1c = HbA1c_pcnt) %>% 
  left_join(cleaned_hw, by = c("id", "date")) %>% 
  mutate(date = as_date(as.numeric(date), origin = "1960-01-01"))

############################ OTHER MEASUREMENTS ################################

convert_flex_ethnicity <- function(eth) {
  if (eth %in% c(1, 2)) {
    return("Not Hispanic or Latino")
  }
  
  if (eth == 3) {
    return("Hispanic or Latino")
  }
  
  if (eth == 4) {
    return("Multiracial/other")
  }
  
  warning("Ethnicity not recognized. Returning `NA`.")
  return(NA)
}

convert_flex_race <- function(race) {
  if (race == 1) {
    return("Black/African American")
  }
  
  if (race == 2) {
    return("American Indian")
  }
  
  if (race == 3) {
    return("Asian")
  }
  
  if (race == 4) {
    return("Hawaiian/Pacific Islander")
  }
  
  if (race == 5) {
    return("White")
  }
  
  if (race == 6) {
    return("Other")
  }
  
  if (race == 7) {
    return("More than one")
  }
  
  warning("Race not recognized. Returning `NA`.")
  return(NA)
}

convert_flex_gender <- function(gender) {
  if (gender == 1) {
    return("F")
  }
  
  if (gender == 2) {
    return("M")
  }
  
  if (gender == 3) {
    return("O")
  }
  
  warning("Gender not recognized. Returning `NA`.")
  return(NA)
}


cleaned_dfs$measurements <-
  dfs$raw_measurements %>% 
  group_by(newID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  transmute(id = newID,
            age = age_yrs,
            years_since_diag = dm_duration,
            ethnicity = map_chr(raceeth_base, convert_flex_ethnicity), 
            race = map_chr(race_youth_base, convert_flex_race),
            gender = map_chr(sex, convert_flex_gender)) %>% 
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

