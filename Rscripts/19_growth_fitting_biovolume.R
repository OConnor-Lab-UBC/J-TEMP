### dealing with biovolume to biomass conversion

## cell dry mass = 0.47(volume)^0.99 (Reynolds 2006)
## cell dry mass =   -0.7550226 + 0.99(log V) (Reynolds 2006)
## strathmann: log C = -0.460 + 0.866(logV)
##  C = 0.109 * V ^ 0.991 (Montagnes et al. 1994, cited in Maranon  2008)
##  logC = -2.216407 + 0.991(logV) (Montagnes et al. 1994, cited in Maranon  2008)
## log C = -0.363 + 0.863 (log V) (Verity et al. 1992)
## from Menden Dueer -1.026 + 1.088 (log V) for chlorophytes
### 0.358(volume)^1.088 


### here the goal is to fit the growth trajectories using biovolume

library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)
library(cowplot)
library(minpack.lm)
library(extrafont)
loadfonts()


tt <- read_csv("data-processed/TT_fit_edit.csv")


### now convert to biomass
tt_mass <- tt %>% 
	mutate(cell_biomass_MD = 0.3584378*(cell_volume)^1.088) %>% 
	mutate(cell_biomass_R = 0.47*(cell_volume)^0.99) %>% 
	mutate(cell_biomass_M = 0.109 *(cell_volume)^0.991) %>% 
	mutate(population_biomass_M = cell_biomass_M * cell_density) %>% 
	mutate(population_biovolume = cell_volume * cell_density)


tt_mass %>% 
	ggplot(aes(y = population_biomass_M, x = days)) + geom_point() + 
	facet_grid(temperature ~ rep)

## find starting biovolume

tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_volume)) %>% View



676.81*1000
fits_many_biovolume <- tt_mass %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(population_biovolume ~ K/(1 + (K/676810 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 100000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 5000000000, r = 200),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))

# get summary info
info_biovolume <- fits_many_biovolume %>%
	unnest(fit %>% map(glance))

# get params
params_biovolume <- fits_many_biovolume %>%
	unnest(fit %>% map(tidy))

new_preds_biovolume <- tt_mass %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

preds_many_fit_biovolume <- fits_many_biovolume %>%
	unnest(fit %>% map(augment, newdata = new_preds))

preds4 <- preds_many_fit_biovolume %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds4, temperature < 33)) +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	geom_point(aes(x = days, y = population_biovolume), data = filter(tt_mass, temperature < 33)) + xlab("Time (days)") + ylab("Population biomass (ug C/ml)")

ggsave("figures/growth_trajectories_biovolume_withCI_32C.pdf", width = 10, height = 10)
# get confidence intervals
CI_biovolume <- fits_many_biovolume %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()


p2 <- params_biovolume %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")

write_csv(p2, "data-processed/K-estimates-biomass-edit.csv")

p2 %>% 
	filter(temperature < 31) %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) +
	geom_smooth(method = "lm", color = "black") + 
	geom_point(size = 2, alpha = 0.5) +
	geom_point(size = 2, shape = 1) +
	labs(y = expression ("Ln carrying capacity \n(population biovolume)"~um^3/ml)) +
	# ylab(("Ln carrying capacity \n(population biovolume um" ~ ^3"/ml)")) +
	xlab("Temperature (1/kT)") +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") 
ggsave("figures/K_biovolume.pdf", width = 4, height = 3.5)

p2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)

# term  estimate std.error statistic      p.value  conf.low conf.high
# 1  (Intercept) 4.7646575 0.7436253  6.407336 4.946987e-06 3.2023587 6.3269564
# 2 inverse_temp 0.2940902 0.0183481 16.028372 4.228887e-12 0.2555422 0.3326381

### result for carrying capacity estimated from trajectories of biovolume. 
p2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% summary()

p2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	ggplot(aes(x = inverse_temp, y = estimate)) + geom_point() +
	scale_x_reverse() + geom_smooth(method = "lm")

# Call:
# 	lm(formula = log(estimate) ~ inverse_temp, data = .)
# 
# Residuals:
# 	Min       1Q   Median       3Q      Max 
# -0.15319 -0.06933  0.01437  0.06063  0.15807 
# 
# Coefficients:
# 	Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.76466    0.74363   6.407 4.95e-06 ***
# 	inverse_temp  0.29409    0.01835  16.028 4.23e-12 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08909 on 18 degrees of freedom
# Multiple R-squared:  0.9345,	Adjusted R-squared:  0.9309 
# F-statistic: 256.9 on 1 and 18 DF,  p-value: 4.229e-12

k_obs_3_biovolume
	