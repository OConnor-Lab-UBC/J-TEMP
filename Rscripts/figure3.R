

library(tidyverse)
library(lubridate)


k_biomass <- read_csv("data-processed/K-estimates-biomass-edit.csv")
cell_sizes <- read_csv("data-processed/cell_sizes.csv")
nitrate <- read_csv("data-processed/nitrate_processed.csv")

nitrate2 <- nitrate %>% 
	mutate(temp = ifelse(temp == 24, 25, temp)) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temp + 273.15))) %>% 
	filter(species == "TT") %>% 
	filter(temp < 32) %>% 
	unite(col = "unique_id", temp, rep, sep = "_", remove =  FALSE)



cell_sizes2 <- cell_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date > "2016-11-28") %>% 
	mutate(cell_biomass_M = 0.109 *(volume_abd)^0.991) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(volume_abd < 2500) %>% 
	group_by(temperature) %>% 
	sample_n(size = 150, replace = FALSE)


cell_sizes2 %>% 
group_by(temperature, rep) %>%
	# ggplot(aes(x = temperature, y = volume_abd, group = date, color = date)) + geom_point() +
	# geom_smooth(method = "lm") + theme_classic()
	summarise(mean_size = mean(cell_biomass_M)) %>% 
	ungroup() %>%
	do(tidy(lm(mean_size ~ temperature, data = .), conf.int = TRUE)) %>% 
	mutate(beta = (estimate/82.2)*100) %>% View

cell_sizes2 %>% 
	filter(temperature == 5) %>%
	summarise(mean_size = mean(cell_biomass_M))  


cell_means <- cell_sizes2 %>% 
	group_by(inverse_temp, rep, temperature) %>%
	# ggplot(aes(x = temperature, y = volume_abd, group = date, color = date)) + geom_point() +
	# geom_smooth(method = "lm") + theme_classic()
	summarise(mean_size = mean(cell_biomass_M)) %>% 
	ungroup()

cell_means %>% 
	lm(mean_size ~ inverse_temp, data = .) %>% summary()


masses <- read_csv("data-processed/cell-masses-dec1.csv")


masses_mean <- masses %>% 
group_by(inverse_temp, temperature) %>% 
	summarise(mean_size = mean(mean_size)) 


masses %>% 
	lm(mean_size ~ inverse_temp, data = .) %>% summary()

masses %>% 
	lm(mean_size ~ temperature, data = .) %>% 
	tidy(., conf.int = TRUE)


cell_size_plot <- masses %>% 
	ggplot(aes(x = inverse_temp, y = mean_size)) +
	geom_jitter(aes(x = inverse_temp, y = cell_biomass_M, group = inverse_temp), data = cell_sizes2, color = "darkgrey", fill = "grey", size = 1.5, width = 0.1, alpha = 0.7) +
	geom_smooth(method = "lm", size =1, color = "black") +
	theme_bw() + geom_point(size = 4, color = "black", alpha = 0.2) +
	geom_point(size = 4, shape = 1) +
	xlab("Temperature (1/kT)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Arial")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab(bquote('Cell size (ug C '*~cell^-1*')')) +
	theme_bw() +
	ggtitle("Temperature (°C)") +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +ylim(0, 175)
	
nitrate_plot <- nitrate2 %>% 
	ggplot(aes(x = inverse_temp, y = nitrate)) +
	geom_smooth(method = "lm", size =1, color = "black") +
	geom_point(size = 4, alpha = 0.2) + 
	geom_point(size = 4, shape = 1) + 
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


k_biomass_plot <- k_biomass %>% 
	filter(temperature < 32) %>% 
	ggplot(aes(x = inverse_temp, y = log(k_biomass))) +
	geom_smooth(method = "lm", size =1, color = "black") +
	geom_point(size = 4, alpha = 0.2) + 
	geom_point(size = 4, shape = 1) + 
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	ylab(bquote('Ln carrying capacity (ug C '*~mL^-1*')')) +
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


k_biomass %>% 
	filter(temperature < 32) %>% 
	lm(log(k_biomass) ~ inverse_temp, data = .) %>% summary()
	
fig3 <- plot_grid(cell_size_plot, k_biomass_plot, nitrate_plot, labels = c("A)", "B)", "C)"), label_fontface = "plain", ncol = 3, nrow = 1, label_x = 0, hjust = 0)
save_plot("figures/k-temp-figure3-edit.pdf", fig3, nrow = 1, ncol = 3, base_height = 3.3, base_width = 3.5)
