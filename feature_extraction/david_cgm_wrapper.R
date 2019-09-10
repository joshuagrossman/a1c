source("lib/cgm_feature_extraction.R")

cgm_df <- NULL
  
cgm_features <- 
  cgm_df %>% 
  make_id_list %>% 
  map(~ c(id = .x$id, make_cgm_feature_df(.x))) %>% 
  reduce(bind_rows)