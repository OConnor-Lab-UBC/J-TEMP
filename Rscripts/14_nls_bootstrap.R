

library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)
library(cowplot)
library(minpack.lm)

sea <- read_csv("data-processed/sea_processed.csv")

TT <- sea %>% 
	filter(species == "TT") %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")


TT_fit <- TT %>% 
filter(cell_density != 36927) %>%
	filter(cell_density != 33992) %>% 
	filter(cell_density != 9279) %>% 
	filter(cell_density != 8924) %>% 
	filter(cell_density != 5045) %>% 
	distinct(cell_density, cell_volume, days, rep, temperature, .keep_all = TRUE)

TT_fit %>% 
	filter(temperature == 38) %>% View

logistic <- function(days, r, K){
	res <- K/(1 + (K/2200 - 1)*exp(-r*days))
	res
}

sub1 <- TT_fit %>% 
	filter(temperature == 5, rep == 1)


fit <- nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
										 data = sub1,
										 iter = 500,
										 start_lower = c(K = 100, r = 0),
										 start_upper = c(K = 10000, r = 1),
										 supp_errors = 'Y',
										 na.action = na.omit,
										lower = c(K = 100, r = 0),
											upper = c(K = 100000, r = 2),
										control = nls.control(maxiter=1000, minFactor=1/204800000))


fit2 <- nlsLM(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
			data= sub1,  start=list(K = 10000, r = 0.1),
			lower = c(K = 100, r = 0),
			upper = c(K = 100000, r = 2),
			control = nls.control(maxiter=1000, minFactor=1/204800000))

class(fit2)
nlsJack(fit2)$jackCI
nlsBoot(fit2)$bootCI
all.equal(fit, fit2)
nlsBoot


fitted(fit)

fits_many <- TT_fit %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 10000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000, r = 2),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))
	
tt_split <- TT_fit %>% 
	split(.$unique_id)



fit_function <- function(df){
	fit <- nlsLM(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
																 data= df,  start=list(K = 10000, r = 0.1),
																 lower = c(K = 100, r = 0),
																 upper = c(K = 100000, r = 2),
																 control = nls.control(maxiter=1000, minFactor=1/204800000))
return(fit)
	}

results <- tt_split %>% 
	map(fit_function)


fit2 <- results[[1]]


boot_function <- function(x){
	out <- x[[1]]
	nb <- nlstools::nlsBoot(x)
	boots <- data.frame(nb$bootCI)
	colnames(boots) <- c("median", "lower", "upper")
	return(boots)
}

results %>% 
	map(boot_function)

boot_many <- group_by(TT_fit, unique_id) %>% 
	# create 200 bootstrap replicates per curve
	do(., boot = modelr::bootstrap(., n = 1000, id = 'boot_num')) %>% 
	# unnest to show bootstrap number, .id
	unnest() %>% 
	# regroup to include the boot_num
	group_by(., unique_id, boot_num) %>% 
	# run the model using map()
	mutate(fit = map(strap, ~ nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
																					data = data.frame(.),
																					iter = 50,
																					start_lower = c(K = 100, r = 0),
																					start_upper = c(K = 10000, r = 1),
																					supp_errors = 'Y',
																					na.action = na.omit,
																					lower = c(K = 100, r = 0),
																					upper = c(K = 50000, r = 2),
																					control = nls.control(maxiter=1000, minFactor=1/204800000))))

save(boot_many, file =  "boot_many.rdata")

boot_many2 <- load(file =  "boot_many.rdata")

new_preds <- TT_fit %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

# get max and min for each curve
max_min <- group_by(TT_fit, unique_id) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create smoother predictions for unbootstrapped models
preds_many_fits <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds))


preds_many <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds))

### do this to take out only the successful fits
preds_id <- boot_many %>%
	unnest(fit %>% map(tidy)) %>% 
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	distinct(uid)

boots_id <- boot_many %>% 
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	distinct(uid)

preds_many_boot <- boot_many %>%
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	filter(uid %in% preds_id$uid) %>% 
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	ungroup() %>% 
	# group by each value of days and get quantiles
	group_by(., unique_id, days) %>%
	summarise(lwr_CI = quantile(.fitted, 0.025),
						upr_CI = quantile(.fitted, 0.975)) %>%
	ungroup() 

write_csv(preds_many_boot, "data-processed/preds_many_boot.csv")

boots_all <- boot_many %>%
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	filter(uid %in% preds_id$uid) %>% 
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	ungroup()

write_csv(boots_all, "data-processed/boots_all.csv")

length(unique(preds_many_boot$unique_id))


preds_boot <- preds_many_boot %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

	# get summary info
info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_many %>%
	unnest(fit %>% map(tidy))



# get confidence intervals
CI <- fits_many %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()



# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))



params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 33) %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
	facet_wrap( ~ term, scales = "free") 

params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)


params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) +
	geom_point(size = 2, alpha = 0.5) +
	geom_point(aes(x = inverse_temp, y = log(K)), data = fits_cool, color = "red") +
	scale_x_reverse() + geom_smooth(method = "lm", color = "black")


# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))


new_preds <- TT_fit %>% 
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(TT_fit, unique_id) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create new predictions
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds))  



preds3 <- preds2 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))


	ggplot(aes(x = days, y = .fitted), data = preds3) + geom_line() + 
	facet_wrap(~ temperature + rep, labeller = labeller(.multi_line = FALSE), ncol = 5) +
	geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = preds_boot, alpha = .2) +
	geom_point(aes(x = days, y = cell_density), data = TT_fit)


	ggplot() +
		geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = preds_boot, alpha = .2) + 
		geom_line(aes(x = days, y = .fitted), data = preds3) +
		facet_wrap(~ temperature + rep, labeller = labeller(.multi_line = FALSE), ncol = 5) +
		geom_point(aes(x = days, y = cell_density), data = TT_fit)
	
	ggsave("figures/growth_trajectories_boot2_withCI.pdf", width = 10, height = 8)
		