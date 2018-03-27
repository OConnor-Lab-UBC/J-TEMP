

library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)
library(cowplot)
library(minpack.lm)
library(soilphysics)
library(extrafont)
loadfonts()
sea1 <- read_csv("data-processed/sea_processed.csv")
sea <- read_csv("data-processed/sea_processed2.csv")

TT_fit <- sea %>% 
	filter(species == "TT") %>% 
	filter(temperature < 38) %>% 
	mutate(cell_density = ifelse(cell_density == 2200, 1000, cell_density)) %>% 
	filter(cell_density != 19767, cell_density != 30185, cell_density != 23949, cell_density != 5638, cell_density != 6505,
				 cell_density != 14164, cell_density != 13597, cell_density != 14438, cell_density != 14650,
				 cell_density != 15049,cell_density != 14530) %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours, start_time) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")



TT_fit %>% 
	filter(temperature == 25, rep == 4) %>% View
	ggplot(aes(x = start_time, y = cell_density)) + geom_point() +
	facet_wrap( ~ rep) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

TT %>% 
	filter(temperature < 32, cell_volume < 2000, days > 34, days < 36) %>% 
	ggplot(aes(x = temperature, y = cell_volume)) + geom_point() +
	facet_wrap( ~ temperature)

TT %>% 
	filter(temperature < 32, cell_volume < 2000, days > 34, days < 35) %>%
	lm(cell_volume ~ temperature, data = .) %>% summary()

TT %>% 
	ggplot(aes(x = days, y = cell_density)) + geom_point() + 
	facet_wrap(temperature ~ rep, ncol =5)


TT_fit <- TT %>% 
	filter(cell_density != 14650) %>% 
	filter(cell_density != 14438) %>% 
	filter(cell_density != 13597) %>% 
	filter(cell_density != 14530) %>% 
	filter(cell_density != 14164) 

# 
# 	filter(temperature == 25, rep ==2) %>% 
# 	filter(days > 20, days < 25) %>% View
# 	ggplot(aes(x = days, y = cell_density, color = factor(days))) + geom_point(size = 4) +
# 	facet_wrap( ~ rep)
	

TT %>% 
	ggplot(aes(x = days, y = cell_density)) + geom_point() +
	facet_wrap(temperature ~ rep, ncol = 5)

TT1 <- sea1 %>% 
	filter(species == "TT") %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")


TT_fit1 <- TT1 %>% 
	filter(temperature < 38) %>% 
	filter(cell_density != 36927) %>%
	filter(cell_density != 33992) %>% 
	filter(cell_density != 9279) %>% 
	filter(cell_density != 8924) %>% 
	filter(cell_density != 5045) %>% 
	distinct(cell_density, cell_volume, days, rep, temperature, .keep_all = TRUE)

TT_fit <- TT %>% 
	filter(temperature < 38) %>% 
filter(cell_density != 36927) %>%
	filter(cell_density != 33992) %>% 
	filter(cell_density != 9279) %>% 
	filter(cell_density != 8924) %>% 
	filter(cell_density != 5045) %>% 
	distinct(cell_density, cell_volume, days, rep, temperature, .keep_all = TRUE)

write_csv(TT_fit, "data-processed/TT_fit.csv")
write_csv(TT_fit, "data-processed/TT_fit2.csv")

TT_fit <- read_csv("data-processed/TT_fit.csv")
TT_fit %>% 
	filter(temperature == 38) %>% View

logistic <- function(days, r, K){
	res <- K/(1 + (K/2200 - 1)*exp(-r*days))
	res
}

sub1 <- TT_fit %>% 
	filter(temperature == 5, rep == 1)


fit <- nls_multstart(cell_density ~ K/(1 + (K/NO - 1)*exp(-r*days)),
										 data = sub1,
										 iter = 500,
										 start_lower = c(K = 100, r = 0, NO = 10),
										 start_upper = c(K = 10000, r = 1, NO = 2000),
										 supp_errors = 'Y',
										 na.action = na.omit,
										 lower = c(K = 100, r = 0, NO = 10),
										 upper = c(K = 50000, r = 2, NO = 2500),
										 control = nls.control(maxiter=1000, minFactor=1/204800000))
summary(fit)

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

TT_25 <- TT_fit %>% 
	filter(temperature == "25") %>% 
	filter(days < 35)





fits_many1 <- TT_fit1 %>% 
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

fits_many <- TT_fit %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(cell_density ~ K/(1 + (K/1000 - 1)*exp(-r*days)),
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

preds_many_boot <- read_csv("data-processed/preds_many_boot.csv")
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


write_csv(params, "data-processed/multstart_params_edit.csv")

params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 33) %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
	facet_wrap( ~ term, scales = "free") 

params1 <- read_csv("data-processed/multstart_params.csv")

params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ungroup() %>% 
	# ggplot(aes(x = inverse_temp, y = estimate)) + geom_point() + geom_smooth(method = "lm") +
	# scale_x_reverse()
	lm(log(estimate) ~ inverse_temp, data = .) %>% summary
	tidy(., conf.int = TRUE)

write_csv(params, "data-processed/multstart_params.csv")
params <- read_csv("data-processed/multstart_params.csv")
p_hot <- params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature == 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")


p_hot %>% 
	filter(estimate < 25000) %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point()

?ggtitle

params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) +
	geom_smooth(method = "lm", color = "black") +
	geom_point(size = 2, alpha = 0.5) +
	geom_point(size = 2, shape = 1) +
	# geom_point(aes(x = inverse_temp, y = log(estimate)), data = p_hot, color = "black", size = 2, alpha = 0.5) +
 ylab("Ln carrying capacity (cells/ml)") +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ggtitle("Temperature (°C)") +
	theme(text = element_text(size=14, family = "Arial")) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figures/figure2_supplement_w_32.pdf", width = 4, height = 3.5)
ggsave("figures/figure2_no_32.pdf", width = 4, height = 3.5)
# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))


new_preds <- TT_fit %>% 
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

# create new predictions

preds1 <- fits_many1 %>%
	unnest(fit %>% map(augment, newdata = new_preds)) 
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds))  

preds3 <- preds2 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

preds1b <- preds1 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))


	ggplot(aes(x = days, y = .fitted), data = preds3) + geom_line() + 
	facet_wrap(~ temperature + rep, labeller = labeller(.multi_line = FALSE), ncol = 5) +
	geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = preds_boot, alpha = .2) +
	geom_point(aes(x = days, y = cell_density), data = TT_fit)


	ggplot() +
		# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
		geom_line(aes(x = days, y = .fitted), data = filter(preds3, temperature < 33)) +
		# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
		facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
		theme(strip.background = element_rect(colour="white", fill="white")) + 
		theme(text = element_text(size=14, family = "Arial")) +
		# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
		geom_point(aes(x = days, y = cell_density), data = filter(TT_fit, temperature < 33)) +
		xlab("Time (days)") + ylab("Population abundance (cells/ml)")
	ggsave("figures/growth_trajectories_boot2_withCI_32C_new.pdf", width = 10, height = 10)	
	ggsave("figures/growth_trajectories_boot2_withCI_32C_old_new.pdf", width = 10, height = 10)	
	ggsave("figures/growth_trajectories_boot2_withCI_32C.pdf", width = 10, height = 10)
		