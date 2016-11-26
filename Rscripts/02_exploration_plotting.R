#### plotting the J-TEMP data

library(tidyverse)
library(plotrix)

jtemp <- read_csv("data-processed/Jtemp_CR_all.csv")

jtemp %>%
	filter(species == "CR", temperature != 35) %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_hours, group = rep, y = total_biovolume, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature) + ggtitle("Chlamydomonas reinhardtii")

jtemp %>% 
	filter(temperature == 8, species == "CR") %>%
	group_by(rep) %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, group = factor(rep), y = cell_density)) +
	# geom_point(size = 4) +
	geom_point() + 
	geom_line() +
	facet_wrap( ~ temperature)



# bring in the seawater species -------------------------------------------

sea_raw <- read_csv("data-processed/YangelJtemp_all.csv")

unique(sea_raw$temperature)

sea <- sea_raw %>% 
	filter(temperature != "18") %>% 
	mutate(temperature = str_replace(temperature, "24", "25")) %>% 
	mutate(temperature = as.numeric(temperature))


# plot the marine species -------------------------------------------------


sea %>%
	filter(temperature != "18") %>% 
	filter(species == "CH") %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_hours, group = rep, y = cell_density, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature)



# join all the data together ----------------------------------------------


all_species <- bind_rows(jtemp, sea)

all_species %>% 
	filter(species == "TT") %>% 
	arrange(desc(cell_density)) %>% View



# plot it all! ------------------------------------------------------------



all_species %>%
	filter(cell_density != 25415) %>% 
	filter(species %in% c("CH", "TT", "CR"), temperature != 35) %>% 
	mutate(month = month(start_time)) %>%
	mutate(day = day(start_time)) %>% 
	mutate(year = 2016) %>% 
	unite(month_day, month, day, year) %>%
	mutate(month_day = mdy(month_day)) %>% 
	group_by(temperature, month_day, species) %>% 
	summarise_each(funs(mean, std.error), total_biovolume) %>% 
	ggplot(aes(x = month_day, y = mean, group = temperature, color = factor(temperature))) +
geom_point(size = 5) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 1) +
	geom_line() + 
	facet_wrap( ~ species) + 
	xlab("date") + ylab("mean cell density") + 
	theme_minimal() +
	theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "grey", size=1))






	