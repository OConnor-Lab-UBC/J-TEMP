

library(tidyverse)
library(stringr)


data <- read_csv("/Users/Joey/Documents/J-TEMP/data-processed/output_rK_TT_cell_abundance.csv")


obs_sim_data <- read_csv("/Users/Joey/Documents/J-TEMP/data-processed/TT_obs_sim_data_cell_abundance.csv")

obs <- obs_sim_data %>% 
	filter(ID == "16_4") %>% 
	filter(type == "observed")

model <- obs_sim_data %>% 
	filter(ID == "16_4") %>% 
	filter(type == "simulated") %>% 
	mutate(type = str_replace(type, "simulated","model"))


obs_sim_data %>% 
	filter(ID == "16_4") %>% 
	mutate(type = str_replace(type, "simulated","model")) %>% 
	ggplot(aes(x = time, y = P, color = type))  + geom_point(size = 2) + 
	scale_color_manual(values = c("orange", "cadetblue")) +
	geom_point(aes(x = time, y = P), data = obs, size = 6, color = "cadetblue") +
	geom_line(aes(x = time, y = P), data = model, size = 2, color = "orange") +
	theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	ylab("population abundance (cells/ml)") + xlab("time (days)") +
	theme(legend.title=element_blank(), legend.text = element_text(size = 18))
ggsave("figures/k-temp-time-series-plot.pdf")



ggplot() +
	geom_point(data = obs_sim_data, aes(x = time, y = P, color = type)) + 
	facet_wrap( ~ ID) + theme_bw() + ylab("population abundance (cells/ml)") + xlab("days")
# ggsave("figures/time_series_fits_TT.png", width = 12, height = 10)


?scale_alpha_manual
