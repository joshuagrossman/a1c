source("lib/feature_extraction/cgm_feature_extraction.R")

create_cgm_features <- function(cgm_df) {
  cgm_df %>% 
  make_id_list %>% 
  map(~ c(id = .x$id, make_cgm_feature_df(.x))) %>% 
  reduce(bind_rows)
}