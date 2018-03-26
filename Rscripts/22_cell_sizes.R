

library(tidyverse)
library(broom)
library(lubridate)

size_data <- read_csv("data-processed/k-temp-size.csv")
size_data_nov28 <- read_csv("data-processed/k-temp-all-cell-sizes-nov28.csv")
size_data_nov14 <- read_csv("data-processed/k-temp-all-cell-sizes-nov14.csv")
size_data_dec1 <- read_csv("data-processed/k-temp-all-cell-sizes-dec1.csv")
size_data_oct31 <- read_csv("data-processed/k-temp-all-cell-sizes-oct31.csv")
size_data_oct28 <- read_csv("data-processed/k-temp-all-cell-sizes-oct28.csv")

all_sizes <- bind_rows(size_data_oct28, size_data_dec1, size_data_nov14, size_data_nov28, size_data_oct31)


size_data %>% 
	lm(cell_volume ~ temperature, data = .) %>% summary()

library(viridis)
all_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(temperature < 32) %>% 
	group_by(temperature, rep, date) %>% 
	summarise(mean_size = mean(volume_abd)) %>% 
	ungroup() %>% 
	ggplot(aes(x = temperature, y = mean_size, group = date, color = date)) + geom_point() +
	geom_smooth(method = "lm") + theme_bw() + 
	facet_wrap( ~ date)
	

all_sizes %>% 
	group_by(temperature, rep, date) %>% 
	summarise(mean_size = mean(volume_abd)) %>% 
	group_by(date) %>% 
	do(tidy(lm(mean_size ~ temperature, data = .), conf.int = TRUE)) %>% View
	tidy(., conf.int = TRUE) %>% View


size_data_all_2 %>% 
	lm(volume_abd ~ temperature, data = .) %>% summary()

size_data_all_2 %>% 
	ggplot(aes(x = temperature, y = volume_abd)) + geom_point() +
	geom_smooth(method = "lm")

size_data_all_2 %>% 
	group_by(temperature, rep) %>% 
	summarise(mean_size = mean(volume_abd)) %>%
	ungroup() %>% 
	ggplot(aes(x = temperature, y = mean_size)) + geom_point() +
	geom_smooth(method = "lm")
	lm(mean_size ~ temperature, data = .) %>% summary()
	
	size_data_all %>% 
		filter(temperature == 5) %>% 
		ungroup() %>% 
		summarise(mean_size = mean(volume_abd)) 
	
	(-16.241 / 863)*100
	
	TT_fit %>% 
		filter(temperature == 5) %>% View
