


### linking K to traits.


params <- read_csv("data-processed/multstart_params.csv")
nitrate <- read_csv("data-processed/nitrate_processed.csv")
final_time_data <- read_csv("data-processed/CH_TT_chla_biovolume_final_time.csv") %>% 
	filter(species == "TT") %>% 
	rename(temp = temperature)

nitrate2 <- nitrate %>% 
	mutate(temp = ifelse(temp == 24, 25, temp)) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temp + 273.15))) %>% 
	filter(species == "TT") %>% 
	unite(col = "unique_id", temp, rep, sep = "_", remove =  FALSE)


np <- left_join(params, nitrate2, by = "unique_id")


np %>% 
	filter(term == "K") %>% 
	filter(temp < 38) %>% 
	ggplot(aes(x = nitrate, y = estimate, color = factor(temp))) + geom_point() + geom_smooth(method = "lm", color = "black")

np %>% 
	filter(term == "K") %>% 
	filter(temp < 38) %>% 
	lm(estimate ~ nitrate, data =.) %>% summary


size_data <- final_time_data %>% 
	mutate(temp = ifelse(temp == 24, 25, temp)) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temp + 273.15))) %>% 
	unite(col = "unique_id", temp, replicate, sep = "_", remove =  FALSE)


sp <- full_join(size_data, np)


sp %>% 
	filter(term == "K") %>% 
	ggplot(aes(x = obs, y = estimate)) + geom_point() +
	facet_wrap( ~ type, scales = "free") + geom_smooth(method = "lm")


sp %>% 
	filter(temp < 31) %>% 
	filter(type == "cell concentration (per ml)", term == "K") %>% 
	mutate(nitrate_per_cell = nitrate/obs) %>% 
	# filter(unique_id != "25_5") %>% 
	ggplot(aes(x = nitrate_per_cell, y = estimate)) + geom_point() + geom_smooth(method = "lm", color = "black") 

sp %>% 
	filter(temp < 33) %>% 
	filter(type == "cell concentration (per ml)", term == "K") %>% 
	mutate(nitrate_use = 10 - nitrate) %>% 
	mutate(nitrate_per_cell = nitrate_use/obs) %>% 
	mutate(nitrate_remaining_per_cell = nitrate/obs) %>% 
	ggplot(aes(x = nitrate_remaining_per_cell, y = estimate)) +  geom_smooth(method = "lm", color = "black") +
	theme(text = element_text(size=14, family = "Arial")) +
	geom_point(size = 2, alpha = 0.5) +
	geom_point(size = 2, shape = 1) +
	ylab("Carrying capacity (cells/ml)") + xlab("Per capita nitrate remaining (uM/cell)")
ggsave("figures/k-v-nitrate.pdf", width = 5, height = 3.5)

sp %>% 
	filter(temp < 33) %>% 
	filter(type == "cell concentration (per ml)", term == "K") %>% 
	mutate(nitrate_use = 10 - nitrate) %>% 
	mutate(nitrate_per_cell = nitrate_use/obs) %>% 
	mutate(nitrate_remaining_per_cell = nitrate/obs) %>% 
	lm(estimate ~ nitrate_remaining_per_cell, data = .) %>% summary

sp %>% 
	filter(temp < 33) %>% 
	filter(type == "cell concentration (per ml)", term == "K") %>% 
	mutate(nitrate_use = 10 - nitrate) %>% 
	mutate(nitrate_per_cell = nitrate_use/obs) %>% 
	mutate(nitrate_remaining_per_cell = nitrate/obs) %>% 
	lm(estimate ~ nitrate, data = .) %>% summary

sp %>% 
	filter(temp < 31) %>% 
	filter(type == "cell concentration (per ml)", term == "K") %>% 
	mutate(nitrate_per_cell = nitrate/obs) %>% 
	lm(estimate ~ nitrate_per_cell, data =.) %>% summary()

sp %>% 
	filter(temp < 32) %>% 
	filter(term == "K", type == 'chla concentration (ug/l)') %>% 
	ggplot(aes(x = temp, y = obs, color = temp)) + geom_point() + geom_smooth(method = "lm") 


sp %>% 
	filter(temp < 33) %>% 
	filter(term == "K", type == 'chla concentration (ug/l)') %>% 
	ggplot(aes(x = obs, y = estimate, color = temp)) + geom_point() + geom_smooth(method = "lm") 
	

nd <- sp %>% 
	# filter(temp < 32) %>% 
	filter(term == "K", type %in% c('chla concentration (ug/l)', "cell concentration (per ml)")) %>% 
	spread(key = type, value = obs) %>%
	ungroup() %>%
	rename(chla = 'chla concentration (ug/l)',
				 cells = 'cell concentration (per ml)') %>% 
	mutate(chl_per = chla/cells)


nd %>% 
	filter(temp < 32) %>% 
	ggplot(aes(x = inverse_temp, y = log(chla))) + geom_point() +
	geom_smooth(method = "lm") + scale_x_reverse() 


nd %>% 
	filter(temp < 32) %>% 
lm(log(chla) ~ inverse_temp, data = .) %>% summary()
	tidy(., conf.int = TRUE)
	
	
	
	sp %>% 
		filter(temp < 33) %>% 
		filter(term == "K", type == 'cell size (um3)') %>% 
		ggplot(aes(x = obs, y = estimate)) + geom_smooth(method = "lm", color = "black") +
		theme(text = element_text(size=14, family = "Arial")) +
		geom_point(size = 2, alpha = 0.5) +
		geom_point(size = 2, shape = 1) +
		xlab(bquote('Cell biovolume ('*um^3*')')) +
		ylab("Carrying capacity (cells/ml)")
	ggsave("figures/k-v-cell-biovolume.pdf", width = 5, height = 3.5)
	
	
	
	sp %>% 
		filter(temp < 33) %>% 
		filter(term == "K", type %in% c('cell size (um3)', "cell concentration (per ml)")) %>% 
		spread(key = type, value = obs) %>%
		ungroup() %>%
		rename(cell_size = 'cell size (um3)',
					 cells = 'cell concentration (per ml)') %>% 
		mutate(nitrate_per_cell = nitrate/cells) %>% 
		ggplot(aes(x = cell_size, y = nitrate_per_cell)) + geom_point() + geom_smooth(method = "lm")
	