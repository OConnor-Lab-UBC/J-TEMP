#### plotting the J-TEMP data

library(tidyverse)
library(plotrix)

jtemp <- read_csv("data-processed/Jtemp_CR_all.csv")

jtemp %>%
	filter(species == "CR", temperature != 35) %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_hours, group = rep, y = cell_density, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature)

jtemp %>% 
	filter(temperature == 8, species == "CR") %>%
	group_by(rep) %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, group = factor(rep), y = cell_density)) +
	# geom_point(size = 4) +
	geom_point() + 
	geom_line() +
	facet_wrap( ~ temperature)
	