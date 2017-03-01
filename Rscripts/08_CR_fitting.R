### Chlamy R fitting


# load packages -----------------------------------------------------------

library(tidyverse)
library(broom)
library(gridExtra)
library(lubridate)
library(plotrix)
library(stringr)
library(nlstools)
library(simecol)
library(purrr)


# read in data ------------------------------------------------------------

jtemp <- read_csv("data-processed/Jtemp_all.csv")

CR <- jtemp %>% 
	filter(species == "CR")

CR %>% 
	filter(total_biovolume < 10^8) %>% 
	filter(time_since_innoc_days < 20) %>% 
	filter(temperature != 35) %>% 
	ggplot(aes(x = time_since_innoc_days, y = total_biovolume)) + geom_point() +
	facet_wrap( ~ temperature + rep, scales = "free")

## initial biovolume

CR %>% 
	filter(time_since_innoc_days < 1) %>% 
	summarise(mean_biovol = mean(total_biovolume))

CR <- CR %>% 
	mutate(total_biovolume = ifelse(time_since_innoc_days < 1, 1996465, total_biovolume))

CR_snipppet <- CR %>% 
	filter(time_since_innoc_days < 3, time_since_innoc_days > 1) %>%
	mutate(rep = ifelse(is.na(rep), 1, rep)) %>% 
	mutate(time_since_innoc_days = 0) %>% 
	mutate(total_biovolume = 1996465) %>% 
	select(time_since_innoc_days, total_biovolume, temperature, rep)

CR_new <- bind_rows(CR, CR_snipppet)

CR_new <- CR_new %>% 
	filter(!is.na(rep))

all_cells_CR <-	CR_new %>%
	unite(Unique_ID, temperature, rep, remove = FALSE) %>%  
	rename(replicate = rep) %>% 
	rename(P = total_biovolume,
				 days = time_since_innoc_days,
				 ID = Unique_ID) %>% 
	arrange(ID, days) %>% 
	filter(!grepl("NA", ID)) %>% 
	filter(days >=0)

 

CR_new %>% 
	filter(total_biovolume < 10^8) %>% 
	filter(time_since_innoc_days < 20) %>% 
	filter(temperature != 35) %>% 
	ggplot(aes(x = time_since_innoc_days, y = total_biovolume)) + geom_point() +
	facet_wrap( ~ temperature + rep, scales = "free")


CR_split <- split(all_cells_CR, f = all_cells_CR$ID) ## split the CR dataframe into little mini dataframes, one for each replicate - temperature combination. We will then use a map function to apply our fitting function to each dataframe individually.



# models ------------------------------------------------------------------

Parameters <- c(r = 1, K = 10 ^ 9) ## initial parameter guesses
CRmodel <- new("odeModel",
							 main = function (time, init, parms) {
							 	with(as.list(c(init, parms)), {
							 		dp <-  r * P * (1 - (P / K))
							 		list(c(dp))
							 	})
							 },
							 parms = Parameters,
							 times = c(from = 0, to = 64, by = 0.1), # the time interval over which the model will be simulated.
							 init = c(P = 1996465), # starting biovolume
							 solver = "lsoda" #lsoda will be called with tolerances of 1e-9. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

controlfit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
															 debuglevel = 0, fn = ssqOdeModel,
															 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
															 control = list(trace = T)
	)
	
	r <- coef(fittedCRmodel)[1]
	K <- coef(fittedCRmodel)[2]
	ID <- data$ID[1]
	output <- data.frame(ID, r, K)
	return(output)
}


# fit models --------------------------------------------------------------

Parameters <- c(r = 0.05, K = 10 ^ 8)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, K = 10 ^ 2)
UpperBound <- c(r = 10, K = 10 ^ 9) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

output_CR_all <- CR_split %>% 
	map_df(controlfit)

output_CR_all %>% 
	separate(ID, into = c("temperature", "rep"), remove = FALSE) %>% 
	filter(K < 10^9) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature > 7, temperature < 38) %>% 
	ggplot(aes(x = temperature, y = log(K))) + geom_point(size = 5, alpha = 0.5) + theme_bw() +
	ylab("ln carrying capacity (K)") + xlab("temperature (C)") + ggtitle("Chlamydomonas rheinhardtii")

ggsave("figures/chlamy_r_k_temp.png")



# Visualize the fits ------------------------------------------------------

plotsinglefit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	# parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. 
	
	# The PORT algorithm is used for the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values.
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
															 debuglevel = 0, fn = ssqOdeModel,
															 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
															 control = list(trace = T)
	)
	
	# To display the fitted results we need to create a new OdeModel object. Here
	# we duplicate CRmodel and then alter it to use our new fitted parameters.
	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)
	
	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=64, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	
	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	observeddata$type <- "observed"
	observeddata <- rename(observeddata, time = obstime)
	simulateddata <- ysim
	simulateddata$type <- "simulated"
	output <- bind_rows(observeddata, simulateddata) ## plop everything together into one data frame
	
	
	# Alternative: Plot the results of our model fitting directly.
	# biol_plot <- ggplot() +
	# 	geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) + # Observed data are points
	# 	geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
	# 	labs(x = "Time (days)", y = "Algal Biovolume")
	
	# Output the results as a ggplot2 object
	# output <- biol_plot
	return(output)
}

obs_sim_data_CR_all <- CR_split %>% 
	map_df(plotsinglefit, .id = "ID") ## don't run this here


obs_sim_data_CR_all %>% 
	separate(ID, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature > 5, temperature < 38) %>% 
	mutate(type = str_replace(type, "simulated", "fit")) %>% 
	ggplot(aes(x = time, y = P, color = type)) + geom_point() + 
	facet_wrap( ~ ID, scales = "free") + theme_bw() + ylab("biovolume (um3/ml)") + xlab("days")
