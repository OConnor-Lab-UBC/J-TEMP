
library(tidyverse)
library(cowplot)
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


nd2 <- sp %>% 
	# filter(temp < 32) %>% 
	filter(term == "K", type %in% c('chla concentration (ug/l)', "total biovolume concentration (um3/ml)")) %>% 
	spread(key = type, value = obs) %>%
	ungroup() %>%
	rename(chla = 'chla concentration (ug/l)',
				 biovol = 'total biovolume concentration (um3/ml)') 

nd2 %>% 
	ggplot(aes(x = chla, y = biovol, color = factor(temp))) + geom_point() + theme_classic() +
	ylab("Total population biovolume (um3/ml)") + xlab("Chla (ug/L)")


nd %>% 
	filter(temp < 32) %>% 
	ggplot(aes(x = inverse_temp, y = log(chla))) + geom_point() +
	geom_smooth(method = "lm") + scale_x_reverse() 


nd %>% 
	filter(temp < 32) %>% 
lm(log(chla) ~ inverse_temp, data = .) %>% 
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
		ggplot(aes(x = cell_size, y = nitrate_per_cell)) + geom_smooth(method = "lm", color = "black") +
		geom_point(size = 2, alpha = 0.5) +
		geom_point(size = 2, shape = 1) +
		xlab(bquote('Cell biovolume ('*um^3*')')) + ylab("Nitrate remaining per cell (uM/cell)")
	ggsave("figures/nitrate-v-cell-biovolume.pdf", width = 5, height = 3.5)
	
	
TT_fit <-	read_csv("data-processed/TT_fit.csv")
library(plotrix)

TT_fit %>% 
	filter(days < 13, temperature < 32) %>% 
	group_by(temperature, rep) %>% 
	summarise_each(funs(mean, std.error), cell_volume) %>% 
	# ungroup() %>% 
	ggplot(aes(x = temperature, y = cell_volume_mean)) + 
	geom_jitter(width = 0.7, size = 2, alpha = 0.7) + geom_smooth(method = "lm") +
	ylab(bquote('Cell biovolume ('*um^3*')')) + xlab("Temperature (°C)")
	

t5 <- TT_fit %>% 
	filter(days < 13, temperature == 5, days > 11) 

t8 <- TT_fit %>% 
	filter(temperature == 8) %>% 
	filter(days < 8, days > 7) 

t16 <- TT_fit %>% 
	filter(temperature == 16) %>%
	filter(days < 8, days > 7)

t25 <- TT_fit %>% 
	filter(temperature == 25) %>% 
	filter(days < 4, days > 3)

TT_all <- bind_rows(t5, t8, t16, t25) 

TT_all %>% 
	filter(temperature > 5) %>% 
	ggplot(aes(x = temperature, y = cell_volume)) + 
	 geom_smooth(method = "lm", color = "black") +
	geom_jitter(width = 0.7, size = 2, alpha = 0.7) +
	ylab(bquote('Cell biovolume ('*um^3*')')) + xlab("Temperature (°C)")

	

TT_all %>% 
	filter(temperature > 5) %>% 
	lm(cell_volume ~ temperature, data = .) %>% summary()

sp %>% 
	filter(temp < 33) %>% 
	filter(term == "K", type %in% c('cell size (um3)', "cell concentration (per ml)")) %>% 
	spread(key = type, value = obs) %>%
	ungroup() %>%
	rename(cell_size = 'cell size (um3)',
				 cells = 'cell concentration (per ml)') %>% 
	lm(cell_size ~ temp, data = .) %>% summary()


TT_fit %>% 
	filter(temperature < 32) %>% 
group_by(temperature, rep) %>% 
	top_n(n = 3, wt = days) %>% 
	summarise(mean_size = mean(cell_volume)) %>% 
	ggplot(aes(x = temperature, y = mean_size)) + geom_point()


TT_fit %>% 
	filter(temperature < 32) %>% 
	group_by(temperature, rep) %>% 
	top_n(n = 4, wt = days) %>% 
	summarise(mean_size = mean(cell_volume)) %>% 
	lm(mean_size ~ temperature, data = .) %>% summary()



# make figure 3 -----------------------------------------------------------

population_biomass <- sp %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>% 
	filter(temp < 32) %>% 
	filter(type == "total biovolume concentration (um3/ml)") %>% 
	mutate(population_biomass = (0.109*(obs)^0.991)/1000)
population_biomass %>% 
	lm(log(population_biomass) ~ inverse_temp, data = .) %>% summary()
population_biomass %>% 
	lm(log(population_biomass) ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)

cell_size <- sp %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>% 
	filter(temp < 32) %>% 
	filter(type == "cell size (um3)") %>% 
	mutate(cell_biomass = 0.109*(obs)^0.991)

cell_size %>% 
	lm(cell_biomass ~ inverse_temp, data = .) %>% summary()
cell_size %>% 
	lm(cell_biomass ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)

cell_biovolume <- sp %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>% 
	filter(temp < 32) %>% 
	filter(type == "cell size (um3)") %>% 
	mutate(cell_biovolume = obs)


nitrate <- sp %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>% 
	filter(temp < 32) %>% 
	distinct(unique_id, .keep_all = TRUE)


nitrate %>% 
	lm(nitrate ~ inverse_temp, data = .) %>% summary()
cell_size %>% 
	lm(cell_biomass ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)


pop_biomass_plot <- population_biomass %>% 
	ggplot(aes(x = inverse_temp, y = log(population_biomass))) +
	geom_smooth(method = "lm", size =1, color = "black") +
	geom_point(size = 2, alpha = 0.2) + 
	geom_point(size = 2, shape = 1) + 
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab(bquote('Ln population biomass (mg C '*~mL^-1*')')) +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	ggtitle("Temperature (°C)")

nitrate_plot <- nitrate %>% 
	ggplot(aes(x = inverse_temp, y = nitrate)) +
	geom_smooth(method = "lm", size =1, color = "black") +
	geom_point(size = 2, alpha = 0.2) + 
	geom_point(size = 2, shape = 1) + 
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab("Nitrate remaining (uM)") +
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
cell_plot <- cell_size %>% 
	ggplot(aes(x = inverse_temp, y = cell_biomass)) +
	geom_smooth(method = "lm", size =1, color = "black") +
	geom_point(size = 2, alpha = 0.2) + 
	geom_point(size = 2, shape = 1) + 
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab(bquote('Cell size (ug C '*~cell^-1*')')) +
	theme_bw() +
ggtitle("Temperature (°C)") +
	theme(text = element_text(size=12, family = "Arial"),
		panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
		plot.title = element_text(hjust = 0.5, size = 12))
fig3 <- plot_grid(nitrate_plot, cell_plot, pop_biomass_plot, labels = c("A)", "B)", "C)"), label_fontface = "plain", ncol = 3, nrow = 1, label_x = 0, hjust = 0)
save_plot("figures/k-temp-figure3.pdf", fig3, nrow = 1, ncol = 3, base_height = 3.3, base_width = 3.5)


