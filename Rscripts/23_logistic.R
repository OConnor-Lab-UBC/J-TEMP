library(tidyverse)
library(broom)
library(modelr)
library(nls.multstart)

sea <- read_csv("data-processed/sea_processed2.csv")

TT_fit <- sea %>% 
	filter(species == "TT") %>% 
	filter(temperature < 34) %>% 
	# filter(cell_density > 1000) %>% 
	# mutate(cell_density = ifelse(cell_density < 2000, 1000, cell_density)) %>% 
	mutate(cell_density = ifelse(cell_density == 2200, 1200, cell_density)) %>% 
	# mutate(cell_density = ifelse(cell_density < 2200 & temperature == 5, 500, cell_density)) %>% 
	filter(cell_density != 19767, cell_density != 30185, cell_density != 23949, cell_density != 5638, cell_density != 6505,
				 cell_density != 14164, cell_density != 13597, cell_density != 14438, cell_density != 14650,
				 cell_density != 15049,cell_density != 14530, cell_density != 5993) %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours, start_time) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")



write_csv(TT_fit, "data-processed/TT_fit_edit.csv")
write_csv(TT_fit, "data-processed/TT_fit_edit-final.csv") ## without 38C

TT_25 <- TT_fit %>% 
	filter(temperature == 25, days < 17)

TT_25 %>% 
	ggplot(aes(x = days, y = cell_density)) + geom_point() +
	facet_wrap( ~ rep)


fits_many <- TT_fit %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(cell_density ~ K/(1 + (K/1200 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 1000,
																								start_lower = c(K = 1000, r = 0),
																								start_upper = c(K = 20000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 50000, r = 2),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))

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


params <- merge(params, CI, by = intersect(names(params), names(CI)))
write_csv(params, "data-processed/params-edit.csv")


params <- read_csv("data-processed/params-edit.csv")
params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
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
	lm(log(estimate) ~ inverse_temp, data = .) %>% summary
tidy(., conf.int = TRUE)


### now try with mass^1/4

all_p <- left_join(params, masses)
write_csv(all_p, "data-processed/K-params-masses.csv")

all_p %>% 
	# separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ungroup() %>% 
	lm(log(estimate*(mean_size^(3/4))) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)

all_p %>% 
	# separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	# filter(estimate < 50000) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ungroup() %>% 
	lm(log(estimate*(68.51255^(3/4))) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)

all_p %>% 
	summarise(mean_all_size = mean(mean_size))


all_p %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate*(mean_size^(3/4))))) +
	geom_smooth(method = "lm", color = "black") +
	geom_point(size = 4, alpha = 0.5) +
	geom_point(size = 4, shape = 1) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab("Ln (carrying capacity (cells/ml) * M^3/4)") +
	theme(text = element_text(size=12, family = "Arial")) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	ggtitle("Temperature (°C)")

ggsave("figures/figure2_mass_exponent.pdf", width = 4, height = 3.5)


params %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) +
	geom_smooth(method = "lm", color = "black") +
	geom_point(size = 4, alpha = 0.5) +
	geom_point(size = 4, shape = 1) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab("Ln carrying capacity (cells/ml)") +
	theme(text = element_text(size=12, family = "Arial")) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	ggtitle("Temperature (°C)")
ggsave("figures/figure2_no_32_edit.pdf", width = 4, height = 3.5)


new_preds <- TT_fit %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds))  

preds3 <- preds2 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))


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



fit_growth <- function(data){
	df <- data
	res <- nls_multstart(cell_density ~ K/(1 + (K/1200 - 1)*exp(-r*days)),
											 data = df,
											 iter = 500,
											 start_lower = c(K = 100, r = 0),
											 start_upper = c(K = 10000, r = 1),
											 supp_errors = 'Y',
											 na.action = na.omit,
											 lower = c(K = 100, r = 0),
											 upper = c(K = 50000, r = 2),
											 control = nls.control(maxiter=1000, minFactor=1/204800000))
	
	expected<-logistic(df$days, coef(res)[2], coef(res)[1])
	rsqr<-1-sum((df$cell_density-expected)^2)/sum((df$cell_density-mean(df$cell_density))^2)
	names(rsqr) <- "rsquared"
	unique_id <- df$unique_id[[1]]
	all <- cbind(rsqr, unique_id)
	return(all)
}

logistic <- function(days, r, K){
	res <- K/(1 + (K/1200 - 1)*exp(-r*days))
	res
}

tt_split <- TT_fit %>% 
	filter(temperature < 32) %>% 
	split(.$unique_id)


all_output1 <- tt_split %>% 
	map_df(fit_growth) 

rsqrs <- all_output1 %>% 
	t(.) %>% 
	data.frame(.) %>% 
	rename(rsq = X1,
				 unique_id = X2)

write_csv(rsqrs, "data-processed/rsq-edit.csv")

View(rsqrs)
rsqrs %>% 
	filter(!grepl("32", unique_id)) %>% 
	mutate(rsq = as.numeric(as.character(rsq))) %>% 
	summarise(mean_r = mean(rsq))

# Fit 32 separately -------------------------------------------------------

TT_32 <- TT_fit %>% 
	filter(temperature == 32) %>% 
	mutate(cell_density = ifelse(cell_density == 2200, 1200, cell_density)) 

fits_32 <- TT_32 %>% 
	group_by(unique_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(cell_density ~ K/(1 + (K/1000 - 1)*exp(-r*days)),
																								data = .x,
																								iter = 500,
																								start_lower = c(K = 100, r = 0),
																								start_upper = c(K = 5000, r = 1),
																								supp_errors = 'N',
																								na.action = na.omit,
																								lower = c(K = 100, r = 0),
																								upper = c(K = 5000, r = 2),
	
																																			control = nls.control(maxiter=1000, minFactor=1/204800000))))




preds32 <- fits_32 %>%
	unnest(fit %>% map(augment, newdata = new_preds))  

preds32 <- preds32 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))


params32 <- fits_32 %>%
	unnest(fit %>% map(tidy))

write_csv(params32, "data-processed/params32.csv")

# get confidence intervals
CI32 <- fits_32 %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()


params32 <- merge(params32, CI32, by = intersect(names(params32), names(CI32)))


fit_growth32 <- function(data){
	df <- data
	res <- nls_multstart(cell_density ~ K/(1 + (K/1200 - 1)*exp(-r*days)),
											 data = df,
											 iter = 500,
											 start_lower = c(K = 100, r = 0),
											 start_upper = c(K = 1000, r = 1),
											 supp_errors = 'Y',
											 na.action = na.omit,
											 lower = c(K = 100, r = 0),
											 upper = c(K = 5000, r = 2),
											 control = nls.control(maxiter=1000, minFactor=1/204800000))
	
	expected<-logistic(df$days, coef(res)[2], coef(res)[1])
	rsqr<-1-sum((df$cell_density-expected)^2)/sum((df$cell_density-mean(df$cell_density))^2)
	names(rsqr) <- "rsquared"
	unique_id <- df$unique_id[[1]]
	all <- cbind(rsqr, unique_id)
	return(all)
}
tt_split32 <- TT_32 %>% 
	split(.$unique_id)


all_output32 <- tt_split32 %>% 
	map_df(fit_growth32) 

rsqrs_32 <- all_output32 %>% 
	t(.) %>% 
	data.frame(.) %>% 
	rename(rsq = X1,
				 unique_id = X2)

ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds32, temperature < 33)) +
	# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
	facet_wrap( ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
	geom_point(aes(x = days, y = cell_density), data = TT_32) +
	xlab("Time (days)") + ylab("Population abundance (cells/ml)")


# Bootstrap 32 ------------------------------------------------------------

boot_32 <- group_by(TT_32, unique_id) %>% 
	# create 200 bootstrap replicates per curve
	do(., boot = modelr::bootstrap(., n = 1000, id = 'boot_num')) %>% 
	# unnest to show bootstrap number, .id
	unnest() %>% 
	# regroup to include the boot_num
	group_by(., unique_id, boot_num) %>% 
	# run the model using map()
	mutate(fit = map(strap, ~ nls_multstart(cell_density ~ K/(1 + (K/1200 - 1)*exp(-r*days)),
																					data = data.frame(.),
																					iter = 50,
																					start_lower = c(K = 100, r = 0),
																					start_upper = c(K = 10000, r = 1),
																					supp_errors = 'Y',
																					na.action = na.omit,
																					lower = c(K = 100, r = 0),
																					upper = c(K = 5000, r = 2),
																					control = nls.control(maxiter=1000, minFactor=1/204800000))))
info32_b <- boot_32 %>%
	unnest(fit %>% map(glance))


preds_id32 <- boot_32 %>%
	unnest(fit %>% map(tidy)) %>% 
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	distinct(uid)

boots_id32 <- boot_32 %>% 
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	distinct(uid)

preds_many_boot32 <- boot_32 %>%
	unite(uid, unique_id, boot_num, remove = FALSE) %>% 
	filter(uid %in% preds_id32$uid) %>% 
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	ungroup() %>% 
	# group by each value of days and get quantiles
	group_by(., unique_id, days) %>%
	summarise(lwr_CI = quantile(.fitted, 0.025),
						upr_CI = quantile(.fitted, 0.975)) %>%
	ungroup() 

write_csv(preds_many_boot32, "data-processed/preds_many_boot_edit32.csv")
preds_many_boot32 <- read_csv("data-processed/preds_many_boot_edit32.csv")

preds_boot32 <- preds_many_boot32 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))

ggplot() +
	geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot32), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds32, temperature < 33)) +
	# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
	geom_point(aes(x = days, y = cell_density), data = TT_32) +
	xlab("Time (days)") + ylab("Population abundance (cells/ml)")

# bootstrap ---------------------------------------------------------------

boot_many <- group_by(TT_fit, unique_id) %>% 
	# create 200 bootstrap replicates per curve
	do(., boot = modelr::bootstrap(., n = 1000, id = 'boot_num')) %>% 
	# unnest to show bootstrap number, .id
	unnest() %>% 
	# regroup to include the boot_num
	group_by(., unique_id, boot_num) %>% 
	# run the model using map()
	mutate(fit = map(strap, ~ nls_multstart(cell_density ~ K/(1 + (K/1200 - 1)*exp(-r*days)),
																					data = data.frame(.),
																					iter = 50,
																					start_lower = c(K = 100, r = 0),
																					start_upper = c(K = 10000, r = 1),
																					supp_errors = 'Y',
																					na.action = na.omit,
																					lower = c(K = 100, r = 0),
																					upper = c(K = 50000, r = 2),
																					control = nls.control(maxiter=1000, minFactor=1/204800000))))

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

write_csv(preds_many_boot, "data-processed/preds_many_boot_edit.csv")
preds_many_boot <- read_csv("data-processed/preds_many_boot_edit.csv")

preds_boot <- preds_many_boot %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))




ggplot() +
	geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = preds_boot32, alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds3, temperature < 32)) +
	geom_line(aes(x = days, y = .fitted), data = filter(preds32, temperature < 33)) +
	# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
	geom_point(aes(x = days, y = cell_density), data = filter(TT_fit, temperature < 33)) +
	xlab("Time (days)") + ylab("Population abundance (cells/ml)")
ggsave("figures/growth_trajectories_withCI_32C_edit_bs.pdf", width = 10, height = 10)



# figure S2 in the supplement ---------------------------------------------


paramscool <- read_csv("data-processed/params-edit.csv") %>% 
	filter(!grepl("32", unique_id)) %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	# filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K")

params32 <- read_csv("data-processed/params32.csv") %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	# filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") 


all_params <- bind_rows(paramscool, params32)

paramscool %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) +
	geom_smooth(method = "lm", color = "black") +
	geom_point(size = 4, alpha = 0.5) +
	geom_point(size = 4, shape = 1) +
	geom_point(aes(x = inverse_temp, y = log(estimate)), data = params32, shape = 1, size = 4) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab("Ln carrying capacity (cells/ml)") +
	theme(text = element_text(size=12, family = "Arial")) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	ggtitle("Temperature (°C)")
ggsave("figures/figure2_with_32_edit.pdf", width = 4, height = 3.5)

log(5000)
