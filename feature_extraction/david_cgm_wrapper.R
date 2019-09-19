source("lib/feature_extraction/cgm_feature_extraction.R")

create_cgm_features <- function(cgm_df) {
  cgm_df %>% 
  make_id_list %>% 
  map(~ c(id = .x$id, make_cgm_feature_df(.x))) %>% 
  # weird hack for init - turn a named vector into a df, then sample 0 rows
  # this ensures that a single patient cgm file will still return a df
  reduce(., bind_rows, .init = sample_n(bind_rows(.[[1]]), 0))
}