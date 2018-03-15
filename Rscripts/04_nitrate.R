#### Nitrate data processing


# load packages -----------------------------------------------------------

library(tidyverse)
library(modelr)
library(gridExtra)
library(cowplot)


# read in data ------------------------------------------------------------


nitrate <- read_csv("data-raw/SO_nitrate_JTEMP.csv")
phosphate <- read_csv("data-raw/SO_phosphate_jtemp.csv")


standards <- nitrate %>% 
	filter(grepl("standard", replicate)) %>% 
	filter(!grepl("COMBO", replicate))
	
	
standards$nitrate <- NA

standards <- standards %>% 
	mutate(nitrate = ifelse(grepl("50uM", replicate), "50", nitrate)) %>% 
	mutate(nitrate = ifelse(grepl("200uM", replicate), "200", nitrate)) %>% 
	mutate(nitrate = ifelse(grepl("500uM", replicate), "500", nitrate)) %>% 
	mutate(nitrate = as.numeric(nitrate))


mod <- lm(nitrate ~ absorbance, data = standards)


ggplot(data = standards, aes(x = nitrate, y = absorbance)) + geom_point() + geom_smooth(method = "lm")


absorbance <- nitrate$absorbance




nitrate_measured <- nitrate %>% 
	add_predictions(mod) %>% 
	rename(nitrate_concentration = pred)


nitrate_processed <- nitrate_measured %>% 
	filter(!grepl("standard", replicate)) %>% 
	separate(replicate, into = c("temperature", "species", "rep"), sep = c(2,4)) %>% 
	mutate(new_temperature = ifelse(grepl("32", temperature), 25, temperature)) %>% 
	mutate(new_temperature = ifelse(grepl("25", temperature), 32, new_temperature)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(new_temperature = as.numeric(new_temperature))


write_csv(nitrate_processed, "data-processed/nitrate_processed_SO.csv")

ggplot(data = nitrate_processed, aes(x = temperature, y = nitrate_concentration)) + geom_point()




# phosphate ---------------------------------------------------------------

standards_phosphate <- phosphate %>% 
	filter(grepl("standard", sample)) %>% 
	filter(!grepl("COMBO", sample)) %>% 
	filter(!grepl("blank", sample)) %>% 
	filter(!is.na(absorbance))

standards_phosphate$phosphate_concentration <- NA

standards_phosphate <- standards_phosphate %>% 
	separate(sample, into = c("sample", "phosphate_standard"), sep = "_") %>% 
	mutate(phosphate_standard = as.numeric(phosphate_standard)) %>% 
	filter(phosphate_standard < 3)


mod1 <- lm(phosphate_standard ~ absorbance, data = standards_phosphate)

df %>% rowwise() %>% mutate(Min = min(A, B, C), Mean = mean(c(A, B, C)))


phosphate_measured <- phosphate %>% 
	rowwise() %>% 
	mutate(absorbance_mean = mean(c(absorbance, absorbance2))) %>% 
	rename(absorbance1 = absorbance) %>% 
	rename(absorbance = absorbance_mean) %>% 
	add_predictions(mod1) %>% 
	rename(phosphate_concentration = pred)

phosphate_samples <- phosphate_measured %>% 
	filter(!grepl("standard", sample)) %>% 
	# filter(!grepl("combo", sample)) %>% 
	separate(sample, into = c("temperature", "species", "rep"), sep = c(2,4)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(species = ifelse(species == "mb", "COMBO", species))
	



write_csv(phosphate_samples, "data-processed/SO_phosphate_concentrations.csv")

ggplot(data = standards_phosphate, aes(x = phosphate_standard, y = absorbance)) + geom_point(size = 3) + geom_smooth(method = "lm")

ggplot(data = phosphate_samples, aes( x = temperature, y = phosphate_concentration)) + geom_point()


### Figure 3 in paper 
library(tidyverse)
library(broom)
nitrate <- read_csv("data-processed/nitrate_processed.csv")
final_time_data <- read_csv("data-processed/CH_TT_chla_biovolume_final_time.csv")

nitrate_stats <- nitrate %>% 
	mutate(temp = ifelse(temp == 24, 25, temp)) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temp + 273.15))) %>% 
	filter(species == "TT") %>% 
	filter(temp < 32) %>% 
	do(tidy(lm(log(nitrate) ~ inverse_temp, data = .), conf.int = TRUE)) 


nitrate_plot <- nitrate %>% 
	mutate(temp = as.numeric(temp)) %>% 
	filter(species == "TT") %>% 
	filter(temp < 32) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temp + 273.15))) %>% 
	ggplot(aes(x = inverse_temp, y = nitrate)) + theme_bw() + xlim(5,25) +
	theme(text = element_text(size = 20)) + ylab("Nitrate (uM N)") +
	# xlab(expression("Temperature (" *degree * "C)")) +
	xlab("") +
	geom_smooth(method = "lm", color = "grey") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + 
	xlab("") +
	ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	geom_point(size = 4, shape = 1, color = "black") +
	geom_point(size = 4, alpha = 0.5)
	annotate("text", x = 16.3, y = 0.5, size = 5, label = "slope = 0.01, 95% CI: -0.09, 0.10")

# ggsave("figures/TT_nitrate_final.png", width = 8, height = 6)


final_time_data %>% 
	filter(species == "TT", type == "cell size (um3)") %>% 
	filter(temperature < 32) %>% 
	do(tidy(lm(obs ~ temperature, data = .), conf.int = TRUE)) 


cell_size_plot <- final_time_data %>% 
	filter(species == "TT", type == "cell size (um3)") %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temperature + 273.15))) %>% 
	ggplot(aes(x = inverse_temp, y = obs)) + theme_bw() +
	theme(text = element_text(size = 20)) + ylab(bquote('Cell size ('*um^3*')')) +
	xlab("") +
	# xlab(expression("Temperature (" *degree * "C)")) +
	geom_smooth(method = "lm", color = "black") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) +
	xlab("") +
	# ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	geom_point(size = 4, shape = 1, color = "black") +
	geom_point(size = 4, alpha = 0.5)
	annotate("text", x = 16.3, y = 1, size = 5, label = "slope = -26.31, 95% CI: -31.64, -20.96")

final_time_data %>% 
	filter(species == "TT", type == "total biovolume concentration (um3/ml)") %>% 
	filter(temperature < 32) %>% 
	mutate(obs = obs/100000) %>% 
	do(tidy(lm(obs ~ temperature, data = .), conf.int = TRUE)) 

final_biovolume <- final_time_data %>% 
	filter(species == "TT", type == "total biovolume concentration (um3/ml)") %>% 
	filter(temperature < 32) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*(temperature + 273.15))) %>% 
	# mutate(obs = obs/100000) %>% 
	ggplot(aes(x = inverse_temp, y = log(obs))) + theme_bw() +
	theme(text = element_text(size = 20)) +  geom_smooth(method = "lm", color = "black") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	xlab("Temperature (1/kT)") + ylab(bquote('ln biovolume ('*um^3~ml^-1*')')) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") +
	# ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	geom_point(size = 4, shape = 1, color = "black") +
	geom_point(size = 4, alpha = 0.5) +
	annotate("text", x = 16.3, y = 1, size = 5, label = "slope = -8.81, 95% CI: -11.10, -6.53")



# p <- grid.arrange(nitrate_plot, cell_size_plot, final_biovolume, nrow =3)
library(cowplot)
figure3 <- plot_grid(nitrate_plot, cell_size_plot, final_biovolume, labels = c("A", "B", "C"), align = "v", ncol = 1, nrow = 3)

save_plot("figures/figure3.png", figure3, nrow = 3, ncol = 1, base_height = 3, base_width = 4)
