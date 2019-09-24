source("lib/feature_extraction/cgm_feature_extraction.R")
library(data.table)
library(lubridate)
library(dplyr)
library(xts)
#########Create function to extarct CGM features#########
create_cgm_features <- function(cgm_df) {
  cgm_df %>% 
    make_id_list %>% 
    map(~ c(id = .x$id, make_cgm_feature_df(.x))) %>% 
    reduce(bind_rows)
}


#########Load CGM data#########
DT1 = fread("~/Downloads/a1c-master/lib/data/OpenAPS/clean/cgm.csv")
DT1$V1 = NULL
DT1$id = paste("Open", DT1$id,  sep="_")
DT1 = unique(DT1)
DT1 = na.omit(DT1)

DT2 = fread("~/Downloads/a1c-master/lib/data/Tidepool/clean/cgm.csv")
DT2$V1 = NULL
DT2$id = paste("Tide", DT2$id, sep="_")
DT2 = unique(DT2)
DT2 = na.omit(DT2)

#########Collapse all data into a list split by patient IDs###########
List1 = split( DT1 , DT1$id )
save(List1, file = "~/Downloads/a1c-master/lib/data/OpenAPS/clean/Open.RDA")

List2 = split( DT2 , DT2$id )
save(List2, file = "~/Downloads/a1c-master/lib/data/Tidepool/clean/Tide.RDA")

#########Concatenate the lists, creating sequential new patient IDs###########
###########Convert datetime into datetime and order by the dt###########
###########Use the dt format to chategorize data for each patient into weeks and years###############
List = append(List1, List2)

List = lapply(1:length(List), function(j) {
####  The use of temp is an artificat of a previous version. It is likely no longer necessary
  temp = List[[j]]
  temp$datetime = ymd_hms(temp$datetime)
#This currently splits by the week of the year. It needs to be updated to split by consecutive two week periods
  temp$id = paste(temp$id, week(temp$datetime), year(temp$datetime), sep = "_")
  temp = temp %>% arrange(temp$datetime)
  return(temp)
})

#### This is where it crashes
Final = lapply(1:length(List), function(j) {
  temp = List[[j]]
  return(create_cgm_features(temp))
})

ff = do.call(rbind, Final)

save(ff, file = "/Users/dscheinker/Downloads/TidePool_and_OpenAPI_List.RDA")

