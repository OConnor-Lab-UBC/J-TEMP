### chla analysis
### last updated by JB jan 5 2017



# load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotrix)


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
	ggplot(aes(x = temperature, y = chla_density, group = species, color = species)) + geom_point() +
	facet_wrap( ~ species)

all_long <- all %>% 
	gather(type, obs, 4:7)
	
all_long %>% 
	filter(species == "TT") %>% 
ggplot(aes(x = temperature, y = obs, color = species, group = species)) + geom_point(size = 4) +
	facet_wrap( ~ type, scales = "free")
