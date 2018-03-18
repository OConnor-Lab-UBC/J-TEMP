

### fit TPC to 2015 data for supplementary materials

library(tidyverse)
library(janitor)


gdata <- read_csv("data-raw/TT_2015.csv") %>% 
	clean_names()


gdata2 <- gdata %>% 
	mutate(days = hours_since_innoc/24)


gdata2 %>% 
	filter(n_treatment == 2) %>% 
	ggplot(aes(x = days, y = particles_per_ml)) + geom_point() +
	facet_wrap(~ temperature)
