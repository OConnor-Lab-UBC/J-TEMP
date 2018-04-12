### new plots for revision
library(ggthemes)
library(broom)
library(tidyverse)

kdata <- read_csv("data-processed/K-params-masses.csv") %>% 
	# separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	mutate(set = "new")

kdata <- read_csv("data-processed/params-edit.csv") %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	mutate(set = "new")


kdata %>% 
	filter(temperature == "5") %>% 
	summarise(mean_size = mean(mean_size))


kdata_hot <- read_csv("data-processed/params32.csv") %>% 
	separate(unique_id, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(term == "K") %>% 
	mutate(set = "new")


kdata_cool <- kdata %>% 
	filter(temperature < 32)



(-15.964/1173)*100
kdata_cool %>% 
	lm(log(estimate*(mean_size^0.75)) ~ inverse_temp, data = .) %>% summary()
kdata_cool %>% 
	lm(log(estimate*(68.51255^0.75)) ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)


kdata_cool %>% 
	lm(log(estimate) ~ inverse_temp, data = .) %>% tidy(., conf.int = TRUE)

kdata_hot <- kdata %>% 
	filter(temperature == 32)

	x <- seq(278.15, 278.15+20, by = 0.01)
	
	tsr_pred<- function(x) {
		y <- (15.5*(81.87 + ((-1.92/100)*81.87)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x)) 
	} 
	
	savage_pred <- function(x) {
		y <- (15.5*(81.87 + ((0/100)*81.87)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
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
	
	plot2a <- ggplot(aes(x = inverse_temp, y = log(estimate*(68.51255^0.75))), data = kdata_cool) + 
		geom_smooth(method = "lm", color = "black") +
		geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
		# geom_smooth(method = "lm", color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed") +
		geom_line(color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed", size =1.5) +
		theme_bw() + geom_point(size = 4, shape = 1, color = "black") +
		geom_point(size = 4, alpha = 0.2) +
		geom_point(data = filter(kdata_hot, log(estimate) < 10), aes(x = inverse_temp, y = log(estimate*(68.51255^0.75))), size = 4, shape = 1) +
		xlab("Temperature (1/kT)") + ylab("Ln (carrying capacity (cells/mL))") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
					panel.background = element_blank(), axis.line = element_line(colour = "black")) +
		theme(text = element_text(size=14, family = "Arial")) +
		scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
		theme(plot.title = element_text(hjust = 0.5, size = 14)) + ylim(12.25, 13.50)
	ggsave("figures/k-temp-figure2_no_32_edit.pdf", width = 5, height = 4)	
	ggsave("figures/k-temp-figure2_with_32_edit.pdf", width = 5, height = 4)	
	
	
	
	tsr_pred_mass <- function(x) {
		y <- 18.5*((86.4 + ((-1.83/100)*86.4)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x)) 
	} 
	
	savage_pred_mass <- function(x) {
		y <- 18.5*((86.4 + ((0/100)*86.4)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
	} 
	
	tsr_predictions_mass <- sapply(x, tsr_pred_mass)
	savage_predictions_mass <- sapply(x, savage_pred_mass)
	
	pred_df_mass <- data.frame(K_tsr = tsr_predictions_mass, K_savage = savage_predictions_mass, temperature_kelvin = x) %>% 
		mutate(inverse_temp = 1/(8.62 * 10^(-5)*temperature_kelvin))
	
	pred_df_mass %>% 
		lm(log(K_savage) ~ inverse_temp, data = .) %>% 
		tidy(., conf.int = TRUE)
	
	
plot2b <- ggplot(aes(x = inverse_temp, y = log(estimate*(mean_size^0.75))), data = kdata_cool) + 
		geom_smooth(method = "lm", color = "black") +
		# geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
		geom_smooth(method = "lm", data = pred_df_mass, aes(x = inverse_temp, y = log(K_savage)), linetype = "dotted", color= "black") +
		theme_bw() + geom_point(size = 4, shape = 1, color = "black") +
		geom_point(size = 4, alpha = 0.2) +
	# geom_abline(slope = 0.33, intercept = -0.399, color = "blue") +
		# geom_point(data = filter(kdata_hot, log(estimate) < 10), aes(x = inverse_temp, y = log(estimate)), size = 4, shape = 1) +
		xlab("Temperature (1/kT)") + ylab("Ln (carrying capacity (cells/mL) * M^3/4)") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
					panel.background = element_blank(), axis.line = element_line(colour = "black")) +
		theme(text = element_text(size=14, family = "Arial")) +
		scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
		theme(plot.title = element_text(hjust = 0.5, size = 14)) + ylim(12.25, 13.50)
	
	
	
fig2 <- plot_grid(plot2a, plot2b, labels = c("A)", "B)"), label_fontface = "plain", ncol = 2, nrow = 1, label_x = 0, hjust = 0)
save_plot("figures/k-temp-figure2-mass-curve.pdf", fig2, nrow = 1, ncol = 2, base_height = 4.2, base_width = 4.5)

	
	pred_df2 %>% 
		lm(log(K_tsr) ~ inverse_temp, data = .) %>% 
		tidy(., conf.int = TRUE)
	
	kdata %>% 
		filter(temperature < 32) %>% 
		lm(log(estimate) ~ inverse_temp, data = .) %>% 
		tidy(., conf.int = TRUE)
	
	pred_df2 %>% 
		lm(log(K_savage) ~ inverse_temp, data = .) %>% 
		tidy(., conf.int = TRUE)
	
	
	
ggsave("figures/k-temp-figure2_w32.pdf", width = 5, height = 4)
ggsave("figures/k-temp-figure2_no_32.pdf", width = 5, height = 4)


library(extrafont)
fonts()


# messing around with equations -------------------------------------------

ggplot(aes(x = inverse_temp, y = log(estimate*(68.51255^0.75))), data = kdata_cool) + 
	# geom_smooth(method = "lm", color = "black") +
	geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
	# geom_smooth(method = "lm", color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed") +
	geom_line(color = "purple", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), size = 1) +
	theme_bw() + 
	# geom_point(size = 4, shape = 1, color = "black") +
	# geom_point(size = 4, alpha = 0.2) +
	xlab("Temperature (1/kT)") + ylab("Ln (carrying capacity (cells/mL))") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Arial")) +
	scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
	theme(plot.title = element_text(hjust = 0.5, size = 14)) 

?scale_x_continuous


h <- function(x) exp(-0.0001*x)

plot(5:25, h(5:25))
