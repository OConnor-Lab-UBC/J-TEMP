#### Nitrate data processing


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data ------------------------------------------------------------


nitrate <- read_csv("data-raw/SO_nitrate_JTEMP.csv")


standards <- nitrate %>% 
	filter(grepl("standard", replicate)) %>% 
	filter(!grepl("COMBO", replicate))
	
	
standards$nitrate <- NA

standards <- standards %>% 
	mutate(nitrate = ifelse(grepl("50uM", replicate), "50", nitrate)) %>% 
	mutate(nitrate = ifelse(grepl("200uM", replicate), "200", nitrate)) %>% 
	mutate(nitrate = ifelse(grepl("500uM", replicate), "500", nitrate)) %>% 
	mutate(nitrate = as.numeric(nitrate))


mod <- lm(nitrate ~ absorbance, data = standards)


ggplot(data = standards, aes(x = nitrate, y = absorbance)) + geom_point() + geom_smooth(method = "lm")


absorbance <- nitrate$absorbance

library(modelr)


nitrate_measured <- nitrate %>% 
	add_predictions(mod) %>% 
	rename(nitrate_concentration = pred)


nitrate_processed <- nitrate_measured %>% 
	filter(!grepl("standard", replicate)) %>% 
	separate(replicate, into = c("temperature", "species", "rep"), sep = c(2,4)) %>% 
	mutate(temperature = as.numeric(temperature))


write_csv(nitrate_processed, "data-processed/nitrate_processed_SO.csv")

ggplot(data = nitrate_processed, aes(x = temperature, y = nitrate_concentration)) + geom_point()



