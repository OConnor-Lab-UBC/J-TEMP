### messing around with the mass dependence stuff
library(tidyverse)
library(broom)

TT <- read_csv("data-processed/output_rK_TT.csv")
size <- read_csv("data-processed/CH_TT_chla_biovolume_final_time.csv")

str(size)
unique(size$type)
size_model <- size %>% 
	filter(type == "cell size (um3)") %>% 
	mutate(inverse_temp = (-1/(.00008617*(temperature+273.15)))) %>%
	filter(species == "TT") %>% 
	do(tidy(lm(obs ~ temperature, data = .), conf.int = TRUE))

TR <- size_model$estimate[size_model$term == "temperature"]

size %>% 
	filter(type == "cell size (um3)") %>% 
	filter(species == "TT") %>% 
	group_by(temperature) %>% 
	summarise(mean_size = mean(obs))

## parameters for calculating changes in size due to increasing temperature.
bR <- -17.34882 # arbitrarily chose a negative number 
M0R <- 1126.092 # base size of resource, at 5C


## effects of temperature on average adult body size ----
# need help writing this function
MR <- M0R + bR*TR #linear function, probably needs to be fine tuned. M0R is baseline adult body size; bR is the slope of the relationship between adult body size and temperature

KMT <- K0*(MR^-3/4)*exp(EP/(k*TR))

K_5 <- 1126.092^-3/4*exp(0.32/((.00008617*(5+273.15))))
K_8 <- 1015.036^-3/4*exp(0.32/((.00008617*(8+273.15))))
K_16 <- 877.248^-3/4*exp(0.32/((.00008617*(16+273.15))))
K_25 <- 575.786^-3/4*exp(0.32/((.00008617*(25+273.15))))

df <- tribble(
	~temperature, ~K,
	5, K_5,
	8, K_8,
	16, K_16,
	25, K_25
)
	
df %>% 
	ggplot(aes(x = temperature, y = K)) + geom_point()

