library(broom)
library(tidyverse)
library(nls.multstart)
library(nlstools)
library(cowplot)
library(minpack.lm)
library(soilphysics)
library(extrafont)
loadfonts()


k38 <- read_csv("data-processed/k-temp-38C.csv") %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")
sea <- read_csv("data-processed/sea_processed2.csv")

TT_38 <- sea %>% 
	filter(species == "TT") %>% 
	filter(temperature == 38) %>% 
	# mutate(cell_density = ifelse(cell_density == 2200, 2200, cell_density)) %>% 
	mutate(cell_density = ifelse(cell_density < 2200, 1000, cell_density)) %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours, start_time) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_") %>% 
	filter(time_since_innoc_hours < 14)


all38 <- bind_rows(k38, TT_38) %>% 
	filter(days < 15)
	mutate(cell_density = ifelse(cell_density < 2200, 1000, cell_density)) %>% 
	mutate(cell_density = ifelse(days > 1, cell_density - 2000, cell_density)) %>% 
	mutate(cell_density = ifelse(cell_density < 0, 0, cell_density))




fits_38 <- all38 %>% 
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
																								upper = c(K = 2000, r = 2),
																								control = nls.control(maxiter=1000, minFactor=1/204800000))))
new_preds <- all38 %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))

preds_38 <- fits_38 %>%
	unnest(fit %>% map(augment, newdata = new_preds))

preds38b <- preds_38 %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature))
info <- fits_38 %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_38 %>%
	unnest(fit %>% map(tidy))



# get confidence intervals
CI <- fits_38 %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
										 	data.frame() %>%
										 	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., unique_id) %>% 
	mutate(., term = c('K', 'r')) %>%
	ungroup()



# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))


ggplot() +
	# geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, x = days), data = filter(preds_boot, temperature < 33), alpha = .3, fill = "grey") + 
	geom_line(aes(x = days, y = .fitted), data = preds38b) +
	# geom_line(aes(x = days, y = .fitted), data = filter(preds1b, temperature < 33), color = "red") +
	facet_wrap( ~ rep, labeller = labeller(.multi_line = FALSE)) +
	theme(strip.background = element_rect(colour="white", fill="white")) + 
	theme(text = element_text(size=14, family = "Arial")) +
	# geom_point(aes(x = days, y = cell_density), data = filter(TT_fit1, temperature < 33), color = "red") +
	geom_point(aes(x = days, y = cell_density), data = all38) +
	xlab("Time (days)") + ylab("Population abundance (cells/ml)") + ylim(0, 30000) + xlim(0, 43)

params %>% 
	ggplot(aes(x = unique_id, y = estimate)) + geom_point() +
	facet_wrap(~ term, scales = "free")


fit_growth <- function(data){
	df <- data
	res <- nls_multstart(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*days)),
											 data = df,
											 iter = 500,
											 start_lower = c(K = 100, r = 0),
											 start_upper = c(K = 10000, r = 1),
											 supp_errors = 'Y',
											 na.action = na.omit,
											 lower = c(K = 100, r = 0),
											 upper = c(K = 2000, r = 2),
											 control = nls.control(maxiter=1000, minFactor=1/204800000))
	
	expected<-logistic(df$days, coef(res)[2], coef(res)[1])
	rsqr<-1-sum((df$cell_density-expected)^2)/sum((df$cell_density-mean(df$cell_density))^2)
	names(rsqr) <- "rsquared"
	unique_id <- df$unique_id[[1]]
	all <- cbind(rsqr, unique_id)
	return(all)
}

logistic <- function(days, r, K){
	res <- K/(1 + (K/2200 - 1)*exp(-r*days))
	res
}



tt_split <- all38 %>% 
	# filter(temperature < 32) %>% 
	split(.$unique_id)


all_output1 <- tt_split %>% 
	map_df(fit_growth) 

rsqrs <- all_output1 %>% 
	t(.) %>% 
	data.frame(.) %>% 
	rename(rsq = X1,
				 unique_id = X2)
