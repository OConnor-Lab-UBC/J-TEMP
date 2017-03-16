library(broom)
library(tidyverse)
library(modelr)


KMT2 <- function(k0, m, EM, EP, k, t) k0*(m^(-3/4))*exp((3*EM - 4*EP)/(4*k*t))

## draw K curve with TSR (in black)
curve(KMT2(k0 = 10000, m = 100, EP = -0.32, EM = -0.01, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

## now add curve for K without TSR, shown in blue
curve(KMT2(k0 = 10000, m = 100, EP = -0.32, EM = 0, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", add = TRUE, col = "blue")

## k0 is the ref K
## m is body mass
## s is the percentage change in body mass with each degree increase in temperature
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis





KMT <- function(k0, m, s, EP, k, t) k0*((m + s*(t-273.15))^(-3/4))*exp(-EP/(k*t))




## draw K curve with TSR (in blue)
curve(KMT(k0 = 10000, m = 500, s = 2.5, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15+5, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "blue", lwd = 3)

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


geom_point(aes(x = temperature_kelvin, y = K), data = k_obs2, size = 3.5, alpha = 0.5) +
### attempt to draw out the predicted curves


KMT <- function(x) 6.5*((1000 + ((-2/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT2 <- function(x) 6.5*((1000 + ((0/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT3 <- function(x) 6.5*((1000 + ((-2.5/100)*1000)*(x-273.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature_kelvin, y = K), data = k_obs2, size = 3.5, alpha = 0.5) +
	stat_function(fun = KMT, color = "black", size = 2) +  stat_function(fun = KMT2, color = "cadetblue", size = 2) +
	stat_function(fun = KMT3, color = "green", size = 2) +
	xlim(273.15, 273.15 + 25) +
	scale_y_continuous(trans = "log", breaks = 5) + theme_bw() +
 xlab("temperature (kelvin)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))



KMT <- function(x) 7.5*((1000 + ((-2/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT2 <- function(x) 7.5*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT3 <- function(x) 7.5*((1000 + ((-2.27/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))


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
ggsave("figures/k-abundance-w-predictions.png")

kpred <- read_csv("data-processed/k-tsr-pred.csv")

kpred2 <- kpred %>% 
	mutate(inverse_temp = 1/(8.62 * 10^(-5)*kelvin))

tsr <- kpred2 %>% 
	filter(treatment == "with_tsr")

no_tsr <- kpred2 %>% 
	filter(treatment == "without_tsr")

mod <- lm(log(k) ~ inverse_temp, data = tsr)
mod_no_tsr <- lm(log(k) ~ inverse_temp, data = no_tsr)

??modelr

kpred3 <- kpred2 %>% 
	add_predictions(mod, var = "pred") 
kpred4 <- kpred2 %>% 
	add_predictions(mod_no_tsr, var = "pred") 

k_obs_3 <- k_obs2 %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 26)

ggplot() + 
	geom_point(aes(x = temperature_kelvin, y = log(K)), data = k_obs_3, size = 4, alpha = 0.5) +
	geom_line(aes(x = kelvin, y = pred-9.5), data = kpred3, color = "black", size = 2) +
	geom_line(aes(x = kelvin, y = pred-9.5), data = kpred4, color = "cadetblue", size = 2) +
	xlim(273.15, 273.15 + 25) + theme_bw() + xlab("temperature (kelvin)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))
	
