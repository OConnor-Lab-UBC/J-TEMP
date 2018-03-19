library(broom)
library(tidyverse)
library(modelr)


## k0 is the ref K
## m is body mass
## s is the percentage change in body mass with each degree increase in temperature
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis




k0 = 100000
m = 1000
s = 0
EP = -0.32
k = 8.62 * 10^(-5)

KMT <- function(k0, m, s, EP, k, t) k0*((m + s*(t-273.15))^(-3/4))*exp(EP/(k*t))

## draw K curve with TSR (in blue)
curve(KMT(k0 = 10000, m = 500, s = 2.5, EP = 0.32, k = 8.62 * 10^(-5), x), from=273.15+5, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "blue", lwd = 3)

RMT <- function(k0, m, s, EP, k, t) k0*((m + s*(t-273.15))^(3/4))*exp(-EP/(k*t))
curve(RMT(k0 = 10, m = 50, s = 2.5, EP = 0.65, k = 8.62 * 10^(-5), x), from=273.15+5, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "blue", lwd = 3)

## draw K curve with TSR (in blue)
curve(KMT(k0 = 10000, m = 500, s = 2.5, EP = 0.65, k = 8.62 * 10^(-5), x), from=273.15+5, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "blue", lwd = 3)




## now add curve for K without TSR, shown in green
curve(KMT(k0 = 10000, m = 500, s = 0, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15+5, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "green", add = TRUE, lwd = 3)



KMT <- function(k0, m, ss, EP, k, t) k0*((m + s*(t - 273.15))^(-3/4))*exp(-EP/(k*t))




## k0 is the ref K
## m is body mass at ref temp
## s is the percentage change in body mass with each degree increase in temperature
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis


KMT <- function(k0, m, s, EP, k, t) k0*((m + ((-s/100)*m)*(t-273.15))^(-3/4))*exp(-EP/(k*t))

## now add curve for K without TSR, shown in blue
curve(KMT(k0 = 100000, m = 1000, s = 0, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15, to=273.15+25, xlab="Temperature (K)", ylab="log(K)", col = "blue", log = "y", lwd = 3)


## draw K curve with TSR (in green)
curve(KMT(k0 = 100000, m = 1000, s = 2.5, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15, to=273.15+25,
			xlab="Temperature (K)", ylab="log(K)", col = "green", log = "y", add = TRUE, lwd = 3)

library(tidyverse)

## now with ggplot!

KMT <- function(k0, m, s, EP, k, x) k0*((m + ((-s/100)*m)*(x-273.15))^(-3/4))*exp(-EP/(k*x))

f <- ggplot(data.frame(x = c(275.15, 275.15+25)), aes(x))

f + stat_function(fun= KMT, k0 = 100000, m = 1000, s = 2.5, EP = -0.32, k = 8.62 * 10^(-5))



## try again!


k0 = 100000
m = 1000
s = 0
EP = -0.32
k = 8.62 * 10^(-5)


KMT <- function(x) 6.5*((1000 + ((-2/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT2 <- function(x) 6.5*((1000 + ((0/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT3 <- function(x) 6.5*((1000 + ((-2.5/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))

p + stat_function(fun = KMT, color = "blue", size = 2) +  stat_function(fun = KMT2, color = "green", size = 2) + xlim(273.15, 273.15 + 25) + scale_y_continuous(trans = "log") + theme_bw() +
	ylab("log K") + xlab("temperature (kelvins)")


k_obs <- read_csv("data-processed/output_rK_TT_cell_abundance.csv")


k_obs2 <- k_obs %>% 
	separate(ID, into = c("temperature", "replicate")) %>%
	mutate(temperature_kelvin = as.numeric(temperature) + 273.15) %>% 
	filter(K < 10^7)



KMT <- function(x) 7.5*((1000 + ((-2/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT2 <- function(x) 7.5*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT3 <- function(x) 7.5*((1000 + ((-2.27/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))

KMT2 <- function(x)17*((1000 + ((0/100)*1000)*(x-278.15))^(-0.87))*exp(0.32/(8.62 * 10^(-5)*x)) ## black
KMT3 <- function(x) 17*((1000 + ((-2.27/100)*1000)*(x-278.15))^(-0.87))*exp(0.32/(8.62 * 10^(-5)*x)) ## grey
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + geom_point(aes(x = temperature_kelvin, y = K), data = k_obs2, size = 0.5, alpha = 0) +
	# stat_function(fun = KMT, color = "grey", size = 2) + 
	stat_function(fun = KMT2, color = "black", size = 3) +
	stat_function(fun = KMT3, color = "grey", size = 3) +
	geom_point(aes(x = temperature_kelvin, y = K), data = k_obs2, size = 6, alpha = 0.5) +
	xlim(273.15 + 5, 273.15 + 25) +
	scale_y_continuous(trans = "log", breaks=seq(8000,30000,10000)) + theme_bw() +
	xlab("temperature (kelvin)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20)) 
ggsave("figures/k-abundance-w-predictions-0.87.png")

kpred <- read_csv("data-processed/k-tsr-pred.csv")

kpreda <- read_csv("data-processed/k-tsr-pred-isometric.csv")
kpreda <- read_csv("data-processed/k-tsr-pred-isometric-0.87.csv")


kpred2 <- kpreda %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*kelvin))

tsr <- kpred2 %>% 
	filter(treatment == "with_tsr")

no_tsr <- kpred2 %>% 
	filter(treatment == "without_tsr")

mod <- lm(log(k) ~ inverse_temp, data = tsr)
mod_no_tsr <- lm(log(k) ~ inverse_temp, data = no_tsr)

summary(mod)
tidy(mod, conf.int = TRUE)
tidy(mod_no_tsr, conf.int = TRUE)

kpred3 <- kpred2 %>% 
	add_predictions(mod, var = "pred") 
kpred4 <- kpred2 %>% 
	add_predictions(mod_no_tsr, var = "pred") 

k_obs_3 <- k_obs2 %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 26) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*temperature_kelvin))

k_obs_4 <- k_obs2 %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*temperature_kelvin))

ggplot() + 
	geom_point(aes(x = temperature_kelvin, y = log(K)), data = k_obs_3, size = 4, alpha = 0.5) +
	geom_line(aes(x = kelvin, y = pred-9.5), data = kpred3, color = "black", size = 2) +
	geom_line(aes(x = kelvin, y = pred-9.5), data = kpred4, color = "cadetblue", size = 2) +
	xlim(273.15, 273.15 + 25) + theme_bw() + xlab("temperature (kelvin)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))



## prediction line
ggplot() +	
ggplot(aes(x = inverse_temp, y = log(K))) + geom_point(size = 6, alpha = 0.5) +geom_smooth(method = "lm", color = "black", size = 3) +
	geom_line(aes(x = inverse_temp, y = pred-9.52), data = kpred4, color = "cadetblue", size = 3) +
	# geom_line(aes(x = inverse_temp, y = pred-9.55), data = kpred3, color = "red", size = 3)+
	scale_x_reverse(limits = c(42, 38.75)) +
	theme_bw() + xlab("temperature (1/kT)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))
ggsave("figures/k-temp-prediction-line.pdf")

## now with data

ggplot(aes(x = inverse_temp, y = log(K)), data = k_obs_3, size = 4, alpha = 0.5) + geom_point(size = 6, alpha = 0.5)+
	geom_smooth(method = "lm", color = "black", size = 3) + 
	geom_line(aes(x = inverse_temp, y = pred-9.52), data = kpred4, color = "cadetblue", size = 3) +
	# geom_line(aes(x = inverse_temp, y = pred-9.55), data = kpred3, color = "red", size = 3)+
	scale_x_reverse(limits = c(42, 38.75)) +
	theme_bw() + xlab("temperature (1/kT)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))
ggsave("figures/k-temp-prediction-line-with-data.pdf")


## now with new prediction (Figure 2 in paper)

ggplot(aes(x = inverse_temp, y = log(K)), data = k_obs_3, size = 4, alpha = 0.5) +
	geom_smooth(method = "lm", color = "#386CB0", size = 0.5) + 
	geom_line(aes(x = inverse_temp, y = pred + 0.025), data = kpred4, color = "black", size = 1, linetype = "dotted") +
	geom_line(aes(x = inverse_temp, y = pred + 0.045), data = kpred3, color = "black", size = 1, linetype = "dashed")+
	geom_smooth(method = "lm", colour="#386CB0", size = 1, fill = "#386CB0") + 
	geom_point(size = 4, alpha = 0.7, colour="#386CB0")+
	geom_point(size = 4, shape = 1, color = "black")+
	# geom_point(aes(x = inverse_temp, y = log(K)), data = k_obs_4, size = 4, alpha = 0.7)+
	scale_x_reverse(limits = c(42, 38.75)) +
	# scale_x_reverse(limits = c(42, 37)) +
	theme_bw() + xlab("Temperature (1/kT)") + ylab("Log carrying capacity (K)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text = element_text(size=14, family = "Helvetica"))
	# annotate("text", x = 41, y = 10.35, label = "prediction incorporating the \n temperature-size rule", size = 5) +
	# annotate("segment", x = 41.75, xend = 41.85, y = 10.34, yend = 10.27, colour="black", size=0.5, arrow = arrow(angle = 15, type = "closed", length = unit(0.03, "npc"))) +
	# annotate("text", x = 40.4, y = 9.4, label = "prediction with temperature \n independent body size", size = 5) +
	# annotate("segment", x = 39.65, xend = 39.4, y = 9.38, yend = 9.38, colour="black", size=0.5, arrow = arrow(angle = 20, type = "closed", length = unit(0.03, "npc"))) 
ggsave("figures/k-temp-prediction-line-with-data-new-pred.pdf")
ggsave("figures/k-temp-prediction-line-with-data-new-pred.png", width = 6, height = 5)


k_obs_biovolume <- read_csv("data-processed/output_rK_TT.csv")

k_obs2_biovolume <- k_obs_biovolume %>% 
	separate(ID, into = c("temperature", "replicate")) %>%
	mutate(temperature_kelvin = as.numeric(temperature) + 273.15) %>% 
	filter(K < 10^9)

k_obs_3_biovolume <- k_obs2_biovolume %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 26) %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*temperature_kelvin))

### now biovolume
ggplot(aes(x = inverse_temp, y = log(K)), data = k_obs_3_biovolume, size = 4, alpha = 0.5) + geom_point(size = 6, alpha = 0.5)+
	geom_smooth(method = "lm", color = "black", size = 3) +
	geom_line(aes(x = inverse_temp, y = pred-2.55), data = kpred4, color = "cadetblue", size = 3) +
	scale_x_reverse(limits = c(42, 38.75)) +
	theme_bw() + xlab("temperature (1/kT)") + ylab("ln population biomass(K)") +
	theme(text = element_text(size=20))
ggsave("figures/k-temp-prediction-line-with-biovolume.pdf")



# now let’s try to actually make predictions ------------------------------

x <- seq(278.15, 278.15+20, by = 0.01)
tsr_pred <- function(x) {
	y <- 5.1*((1000 + ((-2.27/100)*1000)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
} 

savage_pred <- function(x) {
	y <- 4.9*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(0.33/(8.62 * 10^(-5)*x))
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

	ggplot(aes(x = inverse_temp, y = log(K)), data = k_obs_3) + 
		geom_smooth(method = "lm", color = "black") +
	 geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
		geom_smooth(method = "lm", color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed") +
		theme_bw() + geom_point(size = 4, shape = 1, color = "black") +
		geom_point(size = 4, alpha = 0.5) +
		xlab("Temperature (1/kT)") + ylab("ln carrying capacity (cells/mL)") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
					panel.background = element_blank(), axis.line = element_line(colour = "black")) +
		theme(text = element_text(size=14, family = "Helvetica")) +
		scale_x_reverse(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
		theme(plot.title = element_text(hjust = 0.5, size = 14))
	ggsave("figures/k-temp-prediction-line-with-data-dual-axis.png", width = 6, height = 5)


	
## messing around
	

KMT2 <- function(x) 7.5*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(-0.32/(8.62 * 10^(-5)*x))

			
	
		ggplot(aes(x = inverse_temp, y = log(K)), data = k_obs_3) + 
		# geom_smooth(method = "lm", color = "black") +
		# geom_line(aes(x = inverse_temp, y = log(K_savage)), data = pred_df2, linetype = "dotted", size = 1) +
		# geom_smooth(method = "lm", color = "black", data = pred_df2, aes(x = inverse_temp, y = log(K_tsr)), linetype = "dashed") +
		theme_bw() +
		# geom_point(size = 4, shape = 1, color = "black") +
		# geom_point(size = 4, alpha = 0.5) +
		xlab("Temperature (1/kT)") + ylab("ln carrying capacity (cells/mL)") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
					panel.background = element_blank(), axis.line = element_line(colour = "black")) +
		theme(text = element_text(size=14, family = "Helvetica")) +
		scale_x_continuous(sec.axis = sec_axis(~((1/(.*8.62 * 10^(-5)))-273.15))) + xlab("Temperature (1/kT)") + ggtitle("Temperature (°C)") +
		theme(plot.title = element_text(hjust = 0.5, size = 14)) + 
			stat_function(fun = KMT2) + scale_y_log10()
		

