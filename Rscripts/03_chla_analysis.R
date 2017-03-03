### chla analysis
### last updated by JB jan 5 2017



# load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotrix)
library(lubridate)


# read in data ------------------------------------------------------------

chl_raw <- read_csv("data-raw/chla-jtemp.csv")


chl <- chl_raw %>% 
	filter(REPLICATE != "BLANK") %>% 
	filter(REPLICATE != "05tt4") %>% 
	separate(REPLICATE, into = c("temperature", "species", "replicate"), sep = c(2,4)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	rename(rep = replicate) %>% 
	mutate(rep = as.integer(rep))



ggplot(data = chl, aes(x = temperature, y = `chla-calc`, group = species, color = species)) + geom_point()

chl %>% 
	group_by(species, temperature) %>% 
	summarise_each(funs(mean, std.error), `chla-calc`) %>% 
	ggplot(aes(x = temperature, y = mean, group = species, color = species)) + geom_point(size = 4) +
	geom_errorbar(aes(ymin = mean-std.error, ymax = mean + std.error)) + facet_wrap(~ species) +
	scale_y_log10()


# now compare with biovolumes ---------------------------------------------

jtemp <- read_csv("data-processed/sea.csv")

chla <- chl %>% 
	rename(chla = `chla-calc`) %>% 
	select(temperature, rep, species, chla) %>% 
	unite(uniqueid, temperature, rep, species)


dec9 <- jtemp %>% 
	filter(start_time > "2016-12-09") %>% 
	select(temperature, rep, species, cell_density, total_biovolume, cell_volume) %>% 
	unite(uniqueid, temperature, rep, species)

?unite

all <- left_join(chla, dec9, by = "uniqueid") %>% 
	distinct(uniqueid, .keep_all = TRUE) %>% 
	separate(uniqueid, into = c("temperature", "replicate", "species")) %>% 
	mutate(temperature = as.numeric(temperature))

all %>% 
	mutate(chla_density = total_biovolume/chla) %>% 
	mutate(chla_per_cell = cell_density/chla) %>% 
	ggplot(aes(x = temperature, y = chla_density, group = species, color = species)) + geom_point() +
	facet_wrap( ~ species)

all_long <- all %>% 
	mutate(chla_per_biovolume = total_biovolume/chla) %>% 
	mutate(chla_per_cell = cell_density/chla) %>% 
	rename(`cell density (per ml)` = cell_density,
				 `cell size (um3)` = cell_volume,
				 `[chl-a] (ug/l)` = chla,
				 `[biovolume] (um3/ml)` = total_biovolume) %>% 
	gather(type, obs, 4:9)
	
	
all_long %>% 
	filter(species == "TT") %>% 
ggplot(aes(x = temperature, y = obs)) + geom_point(size = 4, alpha = 0.5) +
	facet_wrap( ~ type, scales = "free") + theme_bw() + geom_smooth(method = "lm", color = "black") + ylab("") + xlab("temperature (C)") +
	theme(text = element_text(size=18))
ggsave("figures/TT_cell_physiology.png", width = 8, height = 6) 

all %>% 
	filter(species == "TT") %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
do(tidy(lm(cell_density ~ temperature, data = .), conf.int = TRUE)) %>% View

all %>% 
	filter(species == "TT") %>% 
	lm(cell_volume ~ temperature, data =.) %>% summary()


write_csv(all_long, "data-processed/CH_TT_chla_biovolume_final_time.csv")

### How did cell size change over time?


jtemp %>% 
	filter(species == "TT") %>% 
	filter(cell_volume != 1776.750) %>% 
	mutate(start_time = ifelse(time_since_innoc_days == 0.5, "2016-10-28", start_time)) %>% 
	separate(start_time, into= c("date", "time"), sep = " ") %>% 
	mutate(date = ymd(date)) %>% 
	group_by(temperature, date) %>% 
	summarise_each(funs(mean, std.error), cell_volume) %>% 
	ggplot(aes(x = date, y = mean)) + geom_point(size = 3) +
		geom_errorbar(aes(ymin = mean-std.error, ymax = mean + std.error), width = 0.1) +
	facet_wrap( ~ temperature) + theme_bw() + theme(text = element_text(size = 18)) + ylab("cell size (um3)")
ggsave("figures/TT_cell_size_time_series.png", width = 8, height = 6)

jtemp %>% 
	filter(species == "TT") %>% 
	filter(cell_volume != 1776.750) %>% 
	mutate(start_time = ifelse(time_since_innoc_days == 0.5, "2016-10-28", start_time)) %>% 
	filter(temperature == 38) %>% View
