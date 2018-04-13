

library(tidyverse)
library(broom)
library(lubridate)

size_data <- read_csv("data-processed/k-temp-size.csv")
size_data_nov28 <- read_csv("data-processed/k-temp-all-cell-sizes-nov28.csv")
size_data_nov14 <- read_csv("data-processed/k-temp-all-cell-sizes-nov14.csv")
size_data_dec1 <- read_csv("data-processed/k-temp-all-cell-sizes-dec1.csv")
size_data_oct31 <- read_csv("data-processed/k-temp-all-cell-sizes-oct31.csv")
size_data_oct28 <- read_csv("data-processed/k-temp-all-cell-sizes-oct28.csv")

all_sizes <- bind_rows(size_data_oct28, size_data_dec1, size_data_nov14, size_data_nov28, size_data_oct31)
write_csv(all_sizes, "data-processed/cell_sizes.csv")
all_sizes <- read_csv("data-processed/cell_sizes.csv")

all_sizes$start.time <- ymd("2016-10-28")

all_sizes$day <- interval(all_sizes$start.time, all_sizes$date)
all_sizes$exp_day <- (all_sizes$day/ddays(1)) + 1
all_sizes2 <- all_sizes %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15))))  %>% 
	mutate(exp_day = ifelse(exp_day == 1, "Day 1", exp_day)) %>% 
	mutate(exp_day = ifelse(exp_day == 4, "Day 4", exp_day)) %>% 
	mutate(exp_day = ifelse(exp_day == 18, "Day 18", exp_day)) %>% 
	mutate(exp_day = ifelse(exp_day == 32, "Day 32", exp_day)) %>% 
	mutate(exp_day = ifelse(exp_day == 35, "Day 35", exp_day)) %>% 
	mutate(exp_day = as.factor(exp_day))

all_sizes2$exp_day <- factor(all_sizes2$exp_day, levels = c("Day 1", "Day 4", "Day 18", "Day 32", "Day 35"))

size_data %>% 
	lm(cell_volume ~ temperature, data = .) %>% summary()

library(viridis)
all_sizes2 %>% 
	mutate(date = ymd(date)) %>% 
	filter(temperature < 32) %>% 
	group_by(inverse_temp, rep, exp_day, temperature) %>% 
	summarise(mean_size = mean(volume_abd)) %>% 
	ungroup() %>% 
	mutate(cell_biomass_M = 0.109 *(mean_size)^0.991) %>% 
	ggplot(aes(x = temperature, y = cell_biomass_M))  +
	geom_smooth(method = "lm", color = "black") + theme_classic() + 
	facet_wrap( ~ exp_day) + ylab("Cell biovolume (um3/cell)") + xlab("Temperature (°C)") +
	geom_smooth(method = "lm", size =1, color = "black") +
	theme_bw() + geom_point(size = 4, color = "black", alpha = 0.2) +
	geom_point(size = 4, shape = 1) +
	xlab("Temperature (1/kT)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Arial")) +
	# scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) +
	# scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) +
	xlab("Temperature (°C)") +
	ylab(bquote('Cell size (ug C '*~cell^-1*')')) +
	theme_bw() +
	theme(text = element_text(size=12, family = "Arial"),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "black", size=0.5),
				plot.title = element_text(hjust = 0.5, size = 12)) +
	theme(strip.background = element_rect(colour="white", fill="white")) 
ggsave("figures/cell_size_time.pdf", width = 10, height = 6)
	


all3 <- all_sizes2 %>% 
	mutate(date = ymd(date)) %>% 
	filter(temperature < 32) %>% 
	group_by(inverse_temp, temperature, rep, exp_day) %>% 
	summarise(mean_size = mean(volume_abd)) %>% 
	ungroup() %>% 
	mutate(cell_biomass_M = 0.109 *(mean_size)^0.991) 
	

all3 %>% 
	filter(exp_day == "Day 32") %>% 
	lm(cell_biomass_M ~ temperature, data = .) %>% summary


all_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date == "2016-12-01") %>% 
	mutate(cell_biomass_M = 0.109*(volume_abd)^0.991) %>% 
	mutate(cell_biomass_MD = 0.3584378*(volume_abd)^1.088) %>% 
	group_by(temperature, rep) %>%
	summarise(mean_size = mean(cell_biomass_M)) %>% 
	ungroup() %>%
	do(tidy(lm(mean_size ~ temperature, data = .), conf.int = TRUE)) %>% View


all_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date == "2016-12-01") %>% 
	mutate(cell_biomass_MD = 0.3584378*(volume_abd)^1.088) %>% 
	mutate(cell_biomass_M = 0.109*(volume_abd)^0.991) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	group_by(temperature, rep, inverse_temp) %>%
	# ggplot(aes(x = temperature, y = volume_abd, group = date, color = date)) + geom_point() +
	# geom_smooth(method = "lm") + theme_classic()
	summarise(mean_size = mean(cell_biomass_M)) %>% View
	ungroup() %>%
	lm(mean_size ~ inverse_temp, data = .) %>% summary()

all_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date == "2016-12-01") %>% 
	mutate(cell_biomass_MD = 0.3584378*(volume_abd)^1.088) %>% 
	mutate(cell_biomass_M = 0.109*(volume_abd)^0.991) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	group_by(temperature, rep, inverse_temp) %>%
	# ggplot(aes(x = temperature, y = volume_abd, group = date, color = date)) + geom_point() +
	# geom_smooth(method = "lm") + theme_classic()
	summarise(mean_size = mean(cell_biomass_M)) %>% 
	ungroup() %>%
	lm(mean_size ~ inverse_temp, data = .) %>% summary()


(-15.42491/842.1866)*100
(-18.99854/842.1866)*100## equivalent to a -1.83 decline in cell size per degree C
(-11.85129/842.1866)*100
	
	
masses <- all_sizes %>% 
		mutate(date = ymd(date)) %>% 
		filter(date == "2016-12-01") %>% 
		mutate(cell_biomass_MD = 0.3584378*(volume_abd)^1.088) %>% 
	mutate(cell_biomass_M = 0.109*(volume_abd)^0.991) %>% 
		mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
		group_by(temperature, rep, inverse_temp) %>%
		summarise(mean_size = mean(cell_biomass_M)) %>% 
	unite(col = unique_id, sep = "_", temperature, rep, remove = FALSE)
	
write_csv(masses, "data-processed/cell-masses-dec1.csv")


masses %>% 
	lm(log(mean_size) ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)
masses %>% 
	lm(log(mean_size) ~ inverse_temp, data = .) %>% summary()

masses %>% 
	lm(mean_size ~ temperature, data = .) %>% summary()


masses %>% 
	ggplot(aes(x = inverse_temp, y = log(mean_size))) + geom_point()

masses %>% 
	ggplot(aes(x = temperature, y = mean_size)) + geom_point()



(-1.5716/81.87)
(-1.935209/81.87)
( -1.207894/81.87)
size_function <- function(x) 89.7285 -1.5716*x

size_function(20)

all_sizes %>% 
	mutate(date = ymd(date)) %>% 
	filter(date > "2016-11-29") %>%
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	# group_by(temperature, rep) %>%
	ggplot(aes(x = temperature, y = volume_abd)) + geom_point() +
	geom_smooth(method = "lm") + theme_classic()


(-14.03984/842.209)*100 ## ok on the last day, we're seeing -1.67% decrease in cell size


size_data_all_2 %>% 
	lm(volume_abd ~ temperature, data = .) %>% summary()

size_data_all_2 %>% 
	ggplot(aes(x = temperature, y = volume_abd)) + geom_point() +
	geom_smooth(method = "lm")

size_data_all_2 %>% 
	group_by(temperature, rep) %>% 
	summarise(mean_size = mean(volume_abd)) %>%
	ungroup() %>% 
	ggplot(aes(x = temperature, y = mean_size)) + geom_point() +
	geom_smooth(method = "lm")
	lm(mean_size ~ temperature, data = .) %>% summary()
	
	size_data_all %>% 
		filter(temperature == 5) %>% 
		ungroup() %>% 
		summarise(mean_size = mean(volume_abd)) 
	
	(-16.241 / 863)*100
	
	TT_fit %>% 
		filter(temperature == 5) %>% View
