

library(tidyverse)


k_biomass <- read_csv("data-processed/K-estimates-biomass-edit.csv")
cell_sizes <- read_csv("data-processed/cell_sizes.csv")
nitrate <- read_csv("data-processed/nitrate_processed.csv")





cell_sizes2 <- cell_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date > "2016-11-29") %>% 
	mutate(cell_biomass_M = 0.109 *(volume_abd)^0.991) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(volume_abd < 2500) 

cell_sizes2 %>% 
group_by(temperature, rep) %>%
	# ggplot(aes(x = temperature, y = volume_abd, group = date, color = date)) + geom_point() +
	# geom_smooth(method = "lm") + theme_classic()
	summarise(mean_size = mean(cell_biomass_M)) %>% 
	ungroup() %>%
	do(tidy(lm(mean_size ~ temperature, data = .), conf.int = TRUE)) %>% 
	mutate(beta = (estimate/86.3)*100) %>% View

cell_sizes2 %>% 
	filter(temperature == 5) %>%
	summarise(mean_size = mean(cell_biomass_M))  



cell_sizes2 %>% 
	ggplot(aes(x = inverse_temp, y = cell_biomass_M)) +
	geom_violin(aes(x = inverse_temp, y = cell_biomass_M, group = inverse_temp), data = cell_sizes2) +
	geom_jitter(width = 0.2, size = 2, alpha = 0.2) + 
	geom_smooth(method = "lm", size =1, color = "black") +
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
