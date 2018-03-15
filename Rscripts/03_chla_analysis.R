### chla analysis
### last updated by JB jan 5 2017



# load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotrix)
library(lubridate)
library(broom)

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

all %>% 
	filter(species == "TT") %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temperature + 273.15))) %>% 
	filter(temperature < 32) %>% 
	# ggplot(aes(x = inverse_temp, y = log(chla))) + geom_point()
	do(tidy(lm(log(chla) ~ inverse_temp, data = .), conf.int = TRUE)) 

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
	# filter(temperature < 33) %>% 
	lm(cell_volume ~ temperature, data =.) %>% summary()

all %>% 
	filter(species == "TT") %>%
	ggplot(aes(x = temperature, y = cell_volume)) + geom_point() + geom_smooth(method = "lm")


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


(-17/750)*100

## k0 is the ref K
## m is body mass at ref temp
## s is the percentage change in body mass with each degree increase in temperature
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis


KMT <- function(k0, m, s, EP, k, t) k0*((m + ((-s/100)*m)*(t-273.15))^(-3/4))*exp(-EP/(k*t))

## now add curve for K without TSR, shown in blue
curve(KMT(k0 = 100000, m = 1000, s = 0, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15+0, to=273.15+25, xlab="Temperature (Kelvins)", col = "blue", lwd = 3)


## draw K curve with TSR (in green)
curve(KMT(k0 = 100000, m = 1000, s = 2.26, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15+0, to=273.15+25,
			xlab="Temperature (Kelvins)", ylab="log(K)", col = "green", add = TRUE, lwd = 3)


kpred <- read_csv("data-processed/k-tsr-pred.csv")

kpred %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*kelvin)) %>% 
	group_by(treatment) %>% 
	do(tidy(lm(log(k) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	



### now try to see how k per mass changes with temp


all2 <- all %>% 
	unite(ID, temperature, replicate, sep = "_", remove = FALSE)

?unite


k_obs <- read_csv("data-processed/output_rK_TT_cell_abundance.csv")


k_mass <- left_join(all2, k_obs)


k_mass %>% 
	mutate(mass_st_k = K/(cell_volume^(-3/4))) %>%
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	do(tidy(lm(log(mass_st_k) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View

	
	
