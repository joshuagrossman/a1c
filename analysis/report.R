# report.R
# Generates most stats and plots in the part 1 report.
# Author: Josh Grossman

library(tidyverse)
library(gridExtra)

train_df <- read_csv("data/train.csv")

cor(train_df$mean_bg_full_day, train_df$a1c_value)
cor(train_df$a1c_value, train_df$percent_in_target_range_full_day)
cor(train_df$mean_bg_full_day, 
    train_df$percent_in_target_range_full_day)

male_plot <- 
  ggplot(train_df, aes(x = a1c_value, color = factor(male))) + 
  geom_density() +
  scale_x_continuous(
    name = "HbA1c",
    breaks = 6:11,
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_discrete(
    name = "Sex",
    labels = c("Female", "Male")
  ) +
  scale_y_continuous(
    name = "Density",
    labels = NULL,
    expand = c(0,0)
  ) + 
  theme(legend.position=c(.8,.8))

black_plot <- 
  ggplot(train_df, aes(x = a1c_value, color = factor(black))) + 
  geom_density() +
  scale_x_continuous(
    name = "HbA1c",
    breaks = 6:11,
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_discrete(
    name = "Race",
    labels = c("White", "Black")
  ) +
  scale_y_continuous(
    name = NULL,
    labels = NULL,
    expand = c(0,0)
  ) + 
  theme(legend.position=c(.8,.8))

hispanic_plot <-
  ggplot(train_df, aes(x = a1c_value, color = factor(hispanic))) + 
  geom_density() +
  scale_x_continuous(
    name = "HbA1c",
    breaks = 6:11,
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_discrete(
    name = "Ethnicity",
    labels = c("Not Hispanic", "Hispanic")
  ) +
  scale_y_continuous(
    name = NULL,
    labels = NULL,
    expand = c(0,0)
  ) + 
  theme(legend.position=c(.75,.8))

grid.arrange(male_plot, black_plot, hispanic_plot, nrow = 1) %>% 
  ggsave("demo-density.pdf", ., width = 10, height = 3)

mean(cleaned_df$black)
mean(cleaned_df$hispanic)

ggplot(train_df, aes(x = age, color = dataset)) + geom_density()

cor(train_df$a1c_value, train_df$age)
  