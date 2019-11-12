# Plots a1c and cgm data.

library(GGally)

datasets <- c(
  "RT_CGM_Randomized_Clinical_Trial",
  "CMetforminDataset",
  "FLEX",
  "Protocol_F",
  "Protocol_H"
)

for (d in datasets) {
  feature_df %>% 
    filter(dataset == d) %>% 
    select(a1c_value, race, gender, age, ethnicity, mean_bg_full_day, 
           sd_bg_full_day, percent_in_target_range_full_day, 
           percent_high_full_day) %>% 
    ggpairs(title = d) %>% 
    ggsave(paste0("lib/analysis/", d, ".pdf"), ., height = 20, width = 20)
}

# CMetformin has the fewest a1c measurements by far. 
feature_df %>% 
  group_by(dataset) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = dataset, y = n)) +
  geom_col() + 
  coord_flip()

# CMetformin has the fewest patients as well. Number of a1c measurements
# appears to pretty strongly correlate with number of patients in study.
feature_df %>% 
  group_by(dataset) %>% 
  summarize(n = length(unique(.data$id))) %>% 
  ggplot(aes(x = dataset, y = n)) +
  geom_col() + 
  coord_flip()

# CMetformin and FLEX used only a week of CGM data per a1c measurement.
# Evident from this plot. The other studies have a wider range of available
# CGM data. 
feature_df %>% 
  ggplot(aes(x = cgm_days, fill = dataset)) +
  geom_histogram(position = "dodge", bins = 10)

# Most studies have multiple a1c measurements per patient. 
feature_df %>% 
  ggplot(aes(x = a1c_count, fill = dataset)) +
  geom_histogram(position = "dodge", bins = 10)

# CMetformin and FLEX studies were conducted on individuals with poor diabetes
# control, so it makes sense that mean BGs are higher for those studies.
# Protocol_F falls somewhere in the middle.
feature_df %>% 
  group_by(id, dataset) %>% 
  summarize(mean_bg = mean(mean_bg_full_day, na.rm = T)) %>% 
  ggplot(aes(x = mean_bg, color = dataset)) +
  geom_density()

# CMetformin and FLEX studies were conducted on individuals with poor diabetes
# control, so it makes sense that mean SD BGs are higher for those studies.
# Again, Protocol_F falls somewhere in the middle.
feature_df %>% 
  group_by(id, dataset) %>% 
  summarize(mean_sd_bg = mean(sd_bg_full_day, na.rm = T)) %>% 
  ggplot(aes(x = mean_sd_bg, color = dataset)) +
  geom_density()

# a1c exhibits a similar pattern to mean and sd BG across studies, but noisier.
feature_df %>% 
  group_by(dataset) %>% 
  ggplot(aes(x = a1c_value, color = dataset)) +
  geom_density()

# CMetformin and FLEX are overwhelmingly unhealthy (as expected).
feature_df %>% 
  group_by(dataset, healthy_a1c) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(healthy_prop = n / sum(n)) %>% 
  ggplot(aes(x = healthy_a1c, y = healthy_prop)) +
  geom_col() + 
  coord_flip() +
  facet_wrap(~ dataset)

# CMetformin and FLEX studies were conducted on younger folks.
# Protocol_F is missing ages entirely.
feature_df %>% 
  group_by(dataset) %>% 
  ggplot(aes(x = age, color = dataset)) +
  geom_density()

# Studies other than Protocol_F are almost completely White
# Protocol_F was explcitly designed to measure racial differences
feature_df %>% 
  group_by(dataset, race) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = race, y = prop)) +
  geom_col() + 
  coord_flip() +
  facet_wrap(~ dataset)

# Only CMetformin and FLEX have more than 10% Hispanic particpants.
feature_df %>% 
  group_by(dataset, ethnicity) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = ethnicity, y = prop)) +
  geom_col() + 
  coord_flip() +
  facet_wrap(~ dataset)