#### plotting the J-TEMP data

library(tidyverse)

jtemp <- read_csv("data-processed/Jtemp_all.csv")

jtemp %>%
	group_by(temperature) %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density, color = factor(temperature))) + geom_point()
