#### Nitrate data processing


# load packages -----------------------------------------------------------

library(tidyverse)
library(modelr)


# read in data ------------------------------------------------------------


nitrate <- read_csv("data-raw/SO_nitrate_JTEMP.csv")
phosphate <- read_csv("data-raw/SO_phosphate_jtemp.csv")


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




nitrate_measured <- nitrate %>% 
	add_predictions(mod) %>% 
	rename(nitrate_concentration = pred)


nitrate_processed <- nitrate_measured %>% 
	filter(!grepl("standard", replicate)) %>% 
	separate(replicate, into = c("temperature", "species", "rep"), sep = c(2,4)) %>% 
	mutate(new_temperature = ifelse(grepl("32", temperature), 25, temperature)) %>% 
	mutate(new_temperature = ifelse(grepl("25", temperature), 32, new_temperature)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(new_temperature = as.numeric(new_temperature))


write_csv(nitrate_processed, "data-processed/nitrate_processed_SO.csv")

ggplot(data = nitrate_processed, aes(x = temperature, y = nitrate_concentration)) + geom_point()




# phosphate ---------------------------------------------------------------

standards_phosphate <- phosphate %>% 
	filter(grepl("standard", sample)) %>% 
	filter(!grepl("COMBO", sample)) %>% 
	filter(!grepl("blank", sample)) %>% 
	filter(!is.na(absorbance))

standards_phosphate$phosphate_concentration <- NA

standards_phosphate <- standards_phosphate %>% 
	separate(sample, into = c("sample", "phosphate_standard"), sep = "_") %>% 
	mutate(phosphate_standard = as.numeric(phosphate_standard)) %>% 
	filter(phosphate_standard < 3)


mod1 <- lm(phosphate_standard ~ absorbance, data = standards_phosphate)

df %>% rowwise() %>% mutate(Min = min(A, B, C), Mean = mean(c(A, B, C)))

?add_predictions
phosphate_measured <- phosphate %>% 
	rowwise() %>% 
	mutate(absorbance_mean = mean(c(absorbance, absorbance2))) %>% 
	rename(absorbance1 = absorbance) %>% 
	rename(absorbance = absorbance_mean) %>% 
	add_predictions(mod1) %>% 
	rename(phosphate_concentration = pred)

phosphate_samples <- phosphate_measured %>% 
	filter(!grepl("standard", sample)) %>% 
	# filter(!grepl("combo", sample)) %>% 
	separate(sample, into = c("temperature", "species", "rep"), sep = c(2,4)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(species = ifelse(species == "mb", "COMBO", species))
	



write_csv(phosphate_samples, "data-processed/SO_phosphate_concentrations.csv")

ggplot(data = standards_phosphate, aes(x = phosphate_standard, y = absorbance)) + geom_point(size = 3) + geom_smooth(method = "lm")

ggplot(data = phosphate_samples, aes( x = temperature, y = phosphate_concentration)) + geom_point()





