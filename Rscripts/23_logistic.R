
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
	# ggplot(aes(x = inverse_temp, y = estimate)) + geom_point() + geom_smooth(method = "lm") +
	# scale_x_reverse()
	lm(log(estimate) ~ inverse_temp, data = .) %>% summary
tidy(., conf.int = TRUE)

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


ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = filter(preds3, temperature > 31)) +
	# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
	facet_grid(temperature ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
	geom_point(aes(x = days, y = cell_density), data = filter(TT_fit, temperature >31)) +
	xlab("Time (days)") + ylab("Population abundance (cells/ml)")

fit_growth <- function(data){
	df <- data
	res <- nls_multstart(cell_density ~ K/(1 + (K/1000 - 1)*exp(-r*days)),
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

tt_split <- TT_fit %>% 
	split(.$unique_id)


all_output1 <- tt_split %>% 
	map_df(fit_growth) 

rsqrs <- all_output1 %>% 
	t(.) %>% 
	data.frame(.) %>% 
	rename(rsq = X1,
				 unique_id = X2)
