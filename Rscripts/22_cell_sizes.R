

library(tidyverse)


size_data <- read_csv("data-processed/k-temp-size.csv")
size_data_all <- read_csv("data-processed/k-temp-all-cell-sizes.csv")

size_data %>% 
	lm(cell_volume ~ temperature, data = .) %>% summary()


size_data %>% 
	ggplot(aes(x = temperature, y = cell_volume)) + geom_point() +
	geom_smooth(method = "lm")


size_data_all %>% 
	lm(volume_abd ~ temperature, data = .) %>% summary()

size_data_all %>% 
	ggplot(aes(x = temperature, y = volume_abd)) + geom_point() +
	geom_smooth(method = "lm")

size_data_all %>% 
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
