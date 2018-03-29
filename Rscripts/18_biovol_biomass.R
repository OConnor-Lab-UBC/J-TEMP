
### dealing with biovolume to biomass conversion

## cell dry mass = 0.47(volume)^0.99 (Reynolds 2006)
## cell dry mass =   -0.7550226 + 0.99(log V) (Reynolds 2006)
## strathmann: log C = -0.460 + 0.866(logV)
##  C = 0.109 * V ^ 0.991 (Montagnes et al. 1994, cited in Maranon  2008)
##  logC = -2.216407 + 0.991(logV) (Montagnes et al. 1994, cited in Maranon  2008)
## log C = -0.363 + 0.863 (log V) (Verity et al. 1992)
## from Menden Dueer -1.026 + 1.088 (log V) for chlorophytes
### 0.358(volume)^1.088 


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
	mutate(population_biomass_MD = cell_biomass_MD * cell_density) %>% 
	mutate(population_biomass_R = cell_biomass_R * cell_density)


tt_mass %>% 
	ggplot(aes(y = population_biomass_M, x = days)) + geom_point() + 
	facet_grid(temperature ~ rep)

## find starting biovolume

tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_biomass_MD))


430*1200
fits_many_biomassMD <- tt_mass %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(population_biomass_MD ~ K/(1 + (K/516000 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 100000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000000000, r = 200),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))
tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_biomass_M)) 


tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_biomass_R))

298.5306*1200
fits_many_biomassR <- tt_mass %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(population_biomass_R ~ K/(1 + (K/358236.7 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 100000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000000000, r = 200),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))
tt_mass %>% 
	filter(days < 1) %>% 
	summarise(mean_biomass = mean(cell_biomass_M))
69.5*1200
fits_many_biomassM <- tt_mass %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(population_biomass_M ~ K/(1 + (K/83400 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 100000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000000000, r = 200),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))

# get summary info
info_biomass2 <- fits_many_biomass2 %>%
	unnest(fit %>% map(glance))

info_biomassR <- fits_many_biomassR %>%
	unnest(fit %>% map(glance))

# get params
params_biomassMD <- fits_many_biomassMD %>%
	unnest(fit %>% map(tidy))

params_biomassR <- fits_many_biomassR %>%
	unnest(fit %>% map(tidy))

params_biomassM <- fits_many_biomassM %>%
	unnest(fit %>% map(tidy))

new_preds <- tt_mass %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

preds_many_fits_biomassMD <- fits_many_biomassMD %>%
	unnest(fit %>% map(augment, newdata = new_preds))

preds_many_fits_biomassR <- fits_many_biomassR %>%
	unnest(fit %>% map(augment, newdata = new_preds))

preds3_biomass2 <- preds_many_fits_biomass2 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

preds3_biomassR <- preds_many_fits_biomassR %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds3, temperature < 33)) +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	geom_point(aes(x = days, y = population_biomass_M), data = filter(tt_mass, temperature < 33)) + xlab("Time (days)") + ylab("Population biomass (ug C/ml)")

ggsave("figures/growth_trajectories_biomass_withCI_32C.pdf", width = 10, height = 10)
# get confidence intervals
CI_biomass2 <- fits_many_biomass2 %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()

CI_biomassR <- fits_many_biomassR %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()


p2_biomassM <- params_biomassM %>% ## montagnes
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")

p2_biomassMD <- params_biomassMD %>% ## montagnes
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")

p2_biomassR <- params_biomassR %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")


p2_biomassR %>% 
	filter(temperature < 31) %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 2, alpha = 0.5) +
	geom_point(data = filter(p2_biomass, temperature < 31), aes(x = inverse_temp, y = log(estimate))) +
	geom_smooth(method = "lm", color = "black") + 
	scale_x_reverse()

p2_biomass %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)
# 
# term  estimate  std.error statistic     p.value  conf.low conf.high
# 1  (Intercept) 2.5236749 0.74264964  3.398204 3.20423e-03 0.9634259 4.0839239
# 2 inverse_temp 0.2931953 0.01832403 16.000593 4.35453e-12 0.2546979 0.3316926

p2_biomassMD <- read_csv("data-processed/logistic_parameters_biomass_menden.csv")
p2_biomassM <- read_csv("data-processed/logistic_parameters_biomass_montagnes.csv")
write_csv(p2_biomass2, "data-processed/logistic_parameters_biomass_menden.csv")


p2_biomassM$conversion <- "Montagnes"
p2_biomassM$rep <- as.character(p2_biomassM$rep)
p2_biomassMD$conversion <- "Menden-Deuer"
p2_biomassMD$rep <- as.character(p2_biomassMD$rep)
p2_biomassR$conversion <- "Reynolds"



p2_biomassM %>% 
	filter(temperature < 32) %>% 
	lm(log(estimate) ~ inverse_temp, data =.) %>% 
	tidy(., conf.int = TRUE)

p2_biomassMD %>% 
	filter(temperature < 32) %>% 
	lm(log(estimate) ~ inverse_temp, data =.) %>% 	tidy(., conf.int = TRUE)

p2_biomassR %>% 
	filter(temperature < 32) %>% 
	lm(log(estimate) ~ inverse_temp, data =.) %>% 	tidy(., conf.int = TRUE)

bind_rows(p2_biomassM, p2_biomassMD, p2_biomassR) %>% 
	filter(temperature < 32) %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 2, alpha = 0.5) +
	facet_wrap( ~ conversion, scales = "free") + 
	geom_smooth(method = "lm", color = "black") + 
	scale_x_reverse() +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab(bquote('Ln carrying capacity (ug C '*~mL^-1*')')) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	ggtitle("Temperature (°C)")
		ggsave("figures/K_biomass_comparison.pdf", width = 10, height = 4)

bind_rows(p2_biomass, p2_biomass2, p2_biomassR) %>% 
	filter(temperature < 31) %>% 
	group_by(conversion) %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View

# # Groups:   conversion [2]
# conversion term         estimate std.error statistic  p.value conf.low conf.high
# <chr>      <chr>           <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
# 	1 menden     (Intercept)     3.98     0.754       5.29 5.01e- 5    2.40      5.57 
# 2 menden     inverse_temp    0.303    0.0186     16.3  3.26e-12    0.264     0.342
# 3 montagnes  (Intercept)     2.52     0.743       3.40 3.20e- 3    0.963     4.08 
# 4 montagnes  inverse_temp    0.293    0.0183     16.0  4.35e-12    0.255     0.332



p2_biomass2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)

# term  estimate  std.error statistic      p.value  conf.low conf.high
# 1  (Intercept) 3.9845066 0.75369602  5.286623 5.013708e-05 2.4010500 5.5679632
# 2 inverse_temp 0.3027081 0.01859658 16.277621 3.258348e-12 0.2636381 0.3417781

p2_biomass %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	summary()

# Call:
# 	lm(formula = log(estimate) ~ inverse_temp, data = .)
# 
# Residuals:
# 	Min       1Q   Median       3Q      Max 
# -0.15327 -0.06950  0.01404  0.06053  0.15747 
# 
# Coefficients:
# 	Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.52367    0.74265   3.398   0.0032 ** 
# 	inverse_temp  0.29320    0.01832  16.001 4.35e-12 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08897 on 18 degrees of freedom
# Multiple R-squared:  0.9343,	Adjusted R-squared:  0.9307 
# F-statistic:   256 on 1 and 18 DF,  p-value: 4.355e-12

p2_biomass2 %>% 
	filter(temperature < 31) %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	summary()

# Call:
# 	lm(formula = log(estimate) ~ inverse_temp, data = .)
# 
# Residuals:
# 	Min       1Q   Median       3Q      Max 
# -0.15239 -0.06770  0.01572  0.06150  0.16406 
# 
# Coefficients:	
# 	Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    3.9845     0.7537   5.287 5.01e-05 ***
# 	inverse_temp   0.3027     0.0186  16.278 3.26e-12 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0903 on 18 degrees of freedom
# Multiple R-squared:  0.9364,	Adjusted R-squared:  0.9329 
# F-statistic:   265 on 1 and 18 DF,  p-value: 3.258e-12

write_csv(p2_biomass, "data-processed/logistic_parameters_biomass_montagnes.csv")
write_csv(p2_biomass2, "data-processed/logistic_parameters_biomass_menden.csv")



# plots over time ---------------------------------------------------------

new_preds_biomass <- tt_mass %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

preds_biomassM <- fits_many_biomassM %>%
	unnest(fit %>% map(augment, newdata = new_preds_biomass))

preds_biomassM2 <- preds_biomassM %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds_biomassM2, temperature < 33)) +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	geom_point(aes(x = days, y = population_biomass_M), data = filter(tt_mass, temperature < 33)) + xlab("Time (days)") + ylab("Population biomass (ug C/ml)")
ggsave("figures/growth_trajectories_biomass_32C.pdf", width = 10, height = 10)
