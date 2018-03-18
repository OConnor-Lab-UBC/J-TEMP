

library(tidyverse)
library(minpack.lm)
library(broom)
library(cowplot)


sea <- read_csv("data-processed/sea_processed.csv")


deg5 <- sea %>% 
	filter(temperature == 5, species == "TT") %>% 
	select(rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours))


deg5 %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density)) + geom_point() +
	facet_wrap( ~ rep)
	
	
	fit5 <- nls.lm(cell_density ~ K/(1 + (K/1000 - 1)exp(-r*time_since_innoc_hours)), 
								data = deg5, 
								start = list(K = 10000, r = 0.1),
								lower = c(K = 100, r = 0),
								upper = c(K = 100000, r = 2),
								control = nls.control(maxiter=1024, minFactor=1/204800000))
	
	
	
fit5 <- deg5 %>% 
	group_by(rep) %>% 
 do(tidy(nlsLM(cell_density ~ K/(1 + (K/1000 - 1)*exp(-r*time_since_innoc_hours)),
							 data= deg5,  start=list(K = 10000, r = 0.1),
							 control = nls.control(maxiter=1000, minFactor=1/204800000))))

summary(fit_c)
tidy(fit_c)


TT <- sea %>% 
	filter(species == "TT") %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24)

TT_split <- TT %>% 
	split(.$temperature, .$replicate)


TT_split[[2]]

fit_growth <- function(df){
	res <- try(nlsLM(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
									 data= df,  start=list(K = 10000, r = 0.1),
									 lower = c(K = 100, r = 0),
									 upper = c(K = 100000, r = 2),
									 control = nls.control(maxiter=1000, minFactor=1/204800000)))
	if(class(res)!="try-error"){
		out1 <- tidy(res) %>% 
			select(estimate, term) %>% 
			spread(key = term, value = estimate)
		out2 <- glance(res)
	}
	all <- bind_cols(out1, out2)
	all
}

TT_5 <- TT %>% 
	filter(temperature == 5) %>% 
	split(.$rep)

TT_8 <- TT %>% 
	filter(temperature == 8) %>% 
	split(.$rep)

TT_16 <- TT %>% 
	filter(temperature == 16) %>% 
	filter(cell_density != 36927) %>%
	filter(cell_density != 33992) %>% ## taking out this outlier for now
	split(.$rep)

TT_25 <- TT %>% 
	filter(temperature == 25) %>% 
	split(.$rep)

TT_32 <- TT %>% 
	filter(temperature == 32) %>% 
	split(.$rep)

TT_38 <- TT %>% 
	filter(temperature == 38) %>% 
	split(.$rep)

fits5 <- TT_5 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 5)

fits8 <- TT_8 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 8)

fits16 <- TT_16 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 16)

fits25 <- TT_25 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 25)

fits32 <- TT_32 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 32)

fits38 <- TT_38 %>% 
	map_df(fit_growth, .id = "replicate") %>% 
	mutate(temperature = 38)

all_fits <- bind_rows(fits16, fits25, fits32, fits38, fits5, fits8)


log_fit <- all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	do(tidy(lm(log(K) ~ inverse_temp, data = .), conf.int = TRUE))


all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	lm(log(K) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)


k_fit <- function(x) log_fit$estimate[[1]] + log_fit$estimate[[2]]*x
k_high <- function(x) (log_fit$estimate[[1]] + log_fit$std.error[[1]])  + (log_fit$estimate[[2]] + log_fit$std.error[[2]])*x
k_low <- function(x) (log_fit$estimate[[1]] - log_fit$std.error[[1]])  + (log_fit$estimate[[2]] - log_fit$std.error[[2]])*x


fits_cool <- all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(temperature < 32)
	

fits_hot <- all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(temperature > 31)



fits_cool %>% 
ggplot(aes(x= inverse_temp, y = log(K))) +
	stat_smooth(method = "lm", color = "black") +
	# stat_function(fun = k_fit, color = "red") +
	# stat_function(fun = k_high, color = "red") +
	# stat_function(fun = k_low, color = "red") +
	geom_point(size = 2, alpha = 0.5) +
	# geom_point(size = 2, alpha = 0.5, aes(y = log(K), x = inverse_temp), data = fits_hot) +
	scale_x_reverse() + 
	ylab("Ln(K)") + xlab("Inverse temperature (1/kT)") 
ggsave("figures/ln_K_all_temps.pdf", width = 6, height = 5)


all_fits %>% 
	filter(K < 50000) %>% 
	ggplot(aes(x = temperature, y = r*24)) + geom_point()

TT %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density, group = rep)) + geom_point() +
	facet_wrap(temperature ~ rep)

all_fits2 <- all_fits %>% 
	unite(unique_id, temperature, replicate, remove = FALSE, sep = "_")

reps <- all_fits2 %>% 
	split(.$unique_id)



prediction_function <- function(df){
	pred <- function(x) {
		y <-   df$K/(1 + (df$K/2200 - 1)*exp(-df$r*x))
	}
		x <- seq(0, 43, by = 0.1)
	preds <- sapply(x, pred)
	results <- data.frame(x, preds)
	colnames(results) <- c("time", "abundance")
	return(results)
	
}

predictions <- reps %>% 
	map_df(prediction_function, .id = "unique_id") %>% 
	separate(unique_id, into = c("temperature", "rep")) %>% 
	mutate(temperature = factor(temperature, levels = c("5", "8", "16", "25", "32", "38"))) %>% 
	mutate(part1 = "Temperature = ") %>% 
	mutate(part2 = "Replicate = ") %>% 
	unite(Temperature, part1, temperature, remove = FALSE, sep = "") %>% 
	unite(Rep, part2, rep, remove = FALSE, sep = "") %>% 
	mutate(Temperature = factor(Temperature, 
															levels = c("Temperature = 5", "Temperature = 8", "Temperature = 16", "Temperature = 25", "Temperature = 32", "Temperature = 38"))) 

TT2 <- TT %>% 
	mutate(part1 = "Temperature = ") %>% 
	mutate(part2 = "Replicate = ") %>% 
	unite(Temperature, part1, temperature, remove = FALSE, sep = "") %>% 
	unite(Rep, part2, rep, remove = FALSE, sep = "") %>% 
	mutate(Temperature = factor(Temperature, 
															levels = c("Temperature = 5", "Temperature = 8", "Temperature = 16", "Temperature = 25", "Temperature = 32", "Temperature = 38"))) 
	

	time_series_plot <- TT2 %>% 
		filter(cell_density != 36927) %>%
		filter(cell_density != 33992) %>%
	ggplot(aes(x = days, y = cell_density, group = rep)) + geom_point() +
	geom_line(aes(x = time, y = abundance, group = rep), data = predictions) +
	facet_wrap(Temperature ~ Rep, ncol = 5) + ylab("Population abundance (cells/ml)") + xlab("Time (days)") +
		theme(strip.background = element_rect(colour="white", fill="white")) 
	
	ggsave(time_series_plot, filename = "figures/time_series_facet.pdf", width = 12, height = 10)
	
	

	
	sub1 <- TT2 %>% 
		filter(temperature == 5, rep == 1)
	
	fit_sub1 <- nlsLM(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
				data= sub1,  start=list(K = 10000, r = 0.1),
				lower = c(K = 100, r = 0),
				upper = c(K = 100000, r = 2),
				control = nls.control(maxiter=1000, minFactor=1/204800000))
summary(fit_sub1)
coef(fit_sub1)
str(fit_sub1)
	
	nb <- nlstools::nlsBoot(fit_sub1)
	boots <- data.frame(nb$bootCI)
	colnames(boots) <- c("median", "lower", "upper")

	sub1 %>% 
		ggplot(aes(x = days, y = cell_density)) + geom_point() +
		stat_function(fun = function(x) boots$lower[[1]]/(1 + (boots$lower[[1]]/2200 - 1)*exp(-boots$lower[[2]]*x)), color = "grey")+
		stat_function(fun = function(x) boots$upper[[1]]/(1 + (boots$upper[[1]]/2200 - 1)*exp(-boots$upper[[2]]*x)), color = "grey") +
		stat_function(fun = function(x) boots$median[[1]]/(1 + (boots$median[[1]]/2200 - 1)*exp(-boots$median[[2]]*x)), color = "cadetblue") 
	
	
	
	logistic <- function(days, r, K){
		res <- K/(1 + (K/2200 - 1)*exp(-r*days))
		res
	}
	
	
	expected<-logistic(sub1$days, coef(fit_sub1)[2], coef(fit_sub1)[1])
	rsqr<-1-sum((sub1$cell_density-expected)^2)/sum((sub1$cell_density-mean(sub1$cell_density))^2)
																									
	
	
	### write a function to fit the logistic model, then get the bootstrap CIs and then get the R2s

	## step 1 fit model
	## step 2 calculate Boot CIs
	## step 3 calculate R2
	
	
	