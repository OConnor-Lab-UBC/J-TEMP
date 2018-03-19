
### dealing with biovolume to biomass conversion

## cell dry mass = 0.47(volume)^0.99 (Reynolds 2006)
## cell dry mass =   -0.7550226 + 0.99(log V) (Reynolds 2006)
## strathmann: log C = -0.460 + 0.866(logV)
##  C = 0.0109 * V ^ 0.991 (Montagnes et al. 1994, cited in Maranon  2008)
##  logC = -4.518992 + 0.991(logV) (Montagnes et al. 1994, cited in Maranon  2008)
## log C = -0.363 + 0.863 (log V) (Verity et al. 1992)
## from Menden Dueer -1.026 + 1.088 (log V) for chlorophytes
### 0.358(volume)^1.088 


library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)
library(cowplot)
library(minpack.lm)



tt <- read_csv("data-processed/TT_fit.csv")


### now convert to biomass
tt_mass <- tt %>% 
	mutate(cell_biomass_MD = 0.3584378*(cell_volume)^1.088) %>% 
	mutate(cell_biomass_R = 0.47*(cell_volume)^0.99) %>% 
	mutate(population_biomass_MD = cell_biomass_MD * cell_density)


tt_mass %>% 
	ggplot(aes(y = population_biomass_MD, x = days)) + geom_point() + 
	facet_grid(temperature ~ rep)

## find starting biovolume

tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_biomass_MD)) %>% View


425.1173*2200
fits_many_biomass <- tt_mass %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(population_biomass_MD ~ K/(1 + (K/935258.1 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 100000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000000, r = 200),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))

# get summary info
info <- fits_many_biomass %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_many_biomass %>%
	unnest(fit %>% map(tidy))



# get confidence intervals
CI <- fits_many_biomass %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()


p2 <- params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")


p2 %>% 
	filter(temperature < 31) %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 2, alpha = 0.5) +
	geom_smooth(method = "lm", color = "black") + 
	scale_x_reverse()

p2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE) %>% View
