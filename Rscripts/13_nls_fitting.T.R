

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
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours))

TT_split <- TT %>% 
	split(.$temperature, .$replicate)


TT_split[[2]]

fit_growth <- function(df){
	res <- try(nlsLM(cell_density ~ K/(1 + (K/2200 - 1)*exp(-r*time_since_innoc_hours)),
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

all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ggplot(aes(x= inverse_temp, y = log(K))) + geom_point() +
	scale_x_reverse()


all_fits %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
ungroup() %>% 
	do(tidy(lm(log(K) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View

TT %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density, group = rep)) + geom_point() +
	facet_wrap(~ temperature + rep)


