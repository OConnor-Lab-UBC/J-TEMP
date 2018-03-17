

library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)

sea <- read_csv("data-processed/sea_processed.csv")

TT <- sea %>% 
	filter(species == "TT") %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")


TT_fit <- TT %>% 
filter(cell_density != 36927) %>%
	filter(cell_density != 33992)

logistic <- function(days, r, K){
	res <- K/(1 + (K/2200 - 1)*exp(-r*days))
	res
}

sub1 <- TT2 %>% 
	filter(temperature == 5, rep == 1)


fit <- nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
										 data = sub1,
										 iter = 1000,
										 start_lower = c(K = 100, r = 0),
										 start_upper = c(K = 10000, r = 1),
										 supp_errors = 'N',
										 na.action = na.omit,
										lower = c(K = 100, r = 0),
											upper = c(K = 100000, r = 2),
										control = nls.control(maxiter=1000, minFactor=1/204800000))

fit


fits <- TT_fit %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 1000,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 10000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000, r = 2),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))
	


# get summary info
info <- fits %>%
	unnest(fit %>% map(glance))

# get params
params <- fits %>%
	unnest(fit %>% map(tidy))


params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point() +
	facet_wrap( ~ term, scales = "free")

# get confidence intervals
CI <- fits %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
	unnest(fit %>% map(augment))


new_preds <- TT_fit %>% 
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(TT_fit, unique_id) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create new predictions
preds2 <- fits %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>% 
	merge(., max_min, by = 'curve_id') %>% 
	group_by(., curve_id) %>%
	filter(., days > unique(min_days) & days < unique(max_days)) %>%
	rename(., ln.rate = .fitted) %>%
	ungroup()



preds2 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	ggplot(aes(x = days, y = .fitted)) + geom_line() +
	facet_wrap(~ temperature + rep, labeller = labeller(.multi_line = FALSE)) +
	geom_point(aes(x = days, y = cell_density), data = TT_fit)
