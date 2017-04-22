## nitrate standard curves april 6 2017


library(tidyverse)
library(modelr)


data <- read_csv("data-raw/nitrate-standard-curve-apr06.csv")
data_raw <- read_csv("data-raw/nitrate-uptake-april-07-2017.csv")

### day 2

nitrate_day2 <- data_raw %>% 
	filter(date == "April 7 2017") %>% 
	filter(!grepl("standard", sample)) %>% 
	mutate(absorbance = (absorbance1 + absorbance2 + absorbance3)/3)

standards_day2 <- data_raw %>% 
	filter(date == "April 7 2017") %>% 
	filter(grepl("standard", sample)) %>% 
	separate(sample, into = c("standard", "concentration")) %>% 
	mutate(concentration = as.numeric(concentration)) %>% 
	mutate(absorbance = (absorbance1 + absorbance2 + absorbance3)/3)



mod_day2 <- lm(concentration ~ absorbance, data = standards_day2)

ggplot(data = standards_day2, aes(x = concentration, y = absorbance)) + geom_point() + geom_smooth(method = "lm")

nitrate_measured_day2 <- nitrate_day2 %>% 
	add_predictions(mod_day2) %>% 
	rename(nitrate_concentration = pred) %>% 
	mutate(day = "2")


## day 1


nitrate_day1 <- data_raw %>% 
	filter(date == "April 6 2017") %>% 
	filter(!grepl("standard", sample)) %>% 
	mutate(absorbance = (absorbance1 + absorbance2 + absorbance3)/3)

standards_day1 <- data_raw %>% 
	filter(date == "April 6 2017") %>% 
	filter(grepl("standard", sample)) %>% 
	separate(sample, into = c("standard", "concentration")) %>% 
	mutate(concentration = as.numeric(concentration)) %>% 
	mutate(absorbance = (absorbance1 + absorbance2 + absorbance3)/3)



mod_day1 <- lm(concentration ~ absorbance, data = standards_day1)

ggplot(data = standards_day1, aes(x = concentration, y = absorbance)) + geom_point() + geom_smooth(method = "lm")

nitrate_measured_day1 <- nitrate_day1 %>% 
	add_predictions(mod_day1) %>% 
	rename(nitrate_concentration = pred) %>% 
	mutate(day = "1")



nitrate_all <- bind_rows(nitrate_measured_day1, nitrate_measured_day2)


ggplot(aes(x = day, y = nitrate_concentration, group = sample, color = sample), data = nitrate_all) + geom_point(size = 3) +
	geom_line()

data %>% 
	filter(concentration > 50) %>% 
	filter(time == "afternoon") %>% 
ggplot(aes(x = concentration, y = absorbance)) + geom_point() +
	geom_smooth(method = "lm") +
	theme_bw()

