### new plots for revision
library(ggthemes)
library(broom)


kdata <- read_csv("data-processed/multstart_params_edit.csv") %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	mutate(set = "new")

kdata1 <- read_csv("data-processed/multstart_params.csv") %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	mutate(set = "old")


all_k <- bind_rows(kdata, kdata1)


all_k %>% 
	filter(temperature < 32) %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate), color = set)) + geom_point() +
	geom_smooth(method = "lm") +
	scale_x_reverse()


all_k %>% 
	filter(temperature < 32) %>% 
	group_by(set) %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	tidy(., conf.int = TRUE)


kdata_cool <- kdata %>% 
	filter(temperature < 32)

sea <- read_csv("data-processed/sea_processed2.csv")


sea %>% 
	filter(temperature == 5) %>% 
	filter(time_since_innoc_days > 30) %>% 
	summarise(mean_size = mean(cell_volume))


sea %>% 
	filter(temperature < 32) %>% 
	filter(time_since_innoc_days > 30) %>% 
	group_by(rep, temperature) %>% 
	summarise(mean_size = mean(cell_volume)) %>% 
	lm(mean_size ~ temperature, data = .) %>% summary()


(-15.964/1173)*100



kdata_cool %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)

kdata_cool %>% 
	filter(term == "K") %>% 
	ungroup() %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE) %>% View


kdata_hot <- kdata %>% 
	filter(temperature == 32)

x <- seq(278.15, 278.15+20, by = 0.01)
tsr_pred <- function(x) {
	y <- 4.5*((1000 + ((-1.54/100)*1000)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
} 

savage_pred <- function(x) {
	y <- 4.3*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
} 

tsr_predictions <- sapply(x, tsr_pred)
savage_predictions <- sapply(x, savage_pred)

pred_df <- data.frame(K_tsr = tsr_predictions, K_savage = savage_predictions, temperature_kelvin = x)

pred_df2 <- pred_df %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*temperature_kelvin))

pred_df2 %>% 
	lm(log(K_tsr) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)


pred_df2 %>% 
	lm(log(K_savage) ~ inverse_temp, data = .) %>% 
	tidy(., conf.int = TRUE)


ggplot(aes(x = inverse_temp, y = log(estimate)), data = kdata_cool) + 
	geom_smooth(method = "lm", color = "black") +
	geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
	geom_smooth(method = "lm", color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed") +
	theme_bw() + geom_point(size = 4, shape = 1, color = "black") +
	geom_point(size = 4, alpha = 0.5) +
	geom_point(data = kdata_hot, aes(x = inverse_temp, y = log(estimate)), size = 4, shape = 1) +
	xlab("Temperature (1/kT)") + ylab("Ln carrying capacity (cells/mL)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Arial")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) 
ggsave("figures/k-temp-figure2_w32.pdf", width = 5, height = 4)
ggsave("figures/k-temp-figure2_no_32.pdf", width = 5, height = 4)


library(extrafont)
fonts()
