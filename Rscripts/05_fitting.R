library(tidyverse)
library(simecol)
library(plotrix)
library(gridExtra)
library(stringr)
library(broom)



# bring in the data -------------------------------------------------------

jtemp <- read_csv("data-processed/Jtemp_all.csv")
sea_raw <- read_csv("data-processed/sea.csv")

sea <- sea_raw %>% 
	mutate(total_biovolume = ifelse(is.na(total_biovolume), cell_density*cell_volume, total_biovolume)) %>% 
	filter(cell_volume != 1776.750) %>% 
	filter(cell_volume != 1893.750) %>% 
	filter(cell_volume != 2453.760) %>% 
	filter(cell_volume != 528.400) %>% 
	filter(cell_volume != 581.380) %>% 
	filter(cell_volume != 703.270) %>% 
	filter(cell_volume != 398.070) %>% 
	filter(cell_density > 2100) ## get rid of potentially erroneous measurement points that are really low


write_csv(sea, "data-processed/sea_processed.csv")
	
sea %>% 
	filter(species == "CH", temperature == 38) %>% View

sea %>% 
	filter(species == "CH", temperature == 38) %>% 
	group_by(rep) %>% 
	# filter(time_since_innoc_days < 10) %>% 
	ggplot(aes(x = time_since_innoc_days, y = total_biovolume, color = factor(rep))) + geom_point(size = 2)

sea %>% 
	filter(species == "TT", time_since_innoc_days < 1) %>% 
	group_by(rep, temperature) %>% 
	summarise(mean_biovol = mean(total_biovolume))
	

# set up models -----------------------------------------------------------

CRmodel <- new("odeModel",
							 main = function (time, init, parms) {
							 	with(as.list(c(init, parms)), {
							 		dp <-  r * P * (1 - (P / K))
							 		list(c(dp))
							 	})
							 },
							 parms = Parameters, # Trying Joey's empirically estimated parameters here.
							 times = c(from = 0, to = 43, by = 0.1), # the time interval over which the model will be simulated.
							 init = c(P = 971410.2),
							 solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
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


# plot function -----------------------------------------------------------

plotsinglefit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	# parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. This
	# is fairly standard, although alternatives do exist.
	
	# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values. This may need tweaking.
	
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
	times(plotfittedCRmodel) <- c(from=0, to=43, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	
	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	observeddata$type <- "observed"
	observeddata <- rename(observeddata, time = obstime)
	simulateddata <- ysim
	simulateddata$type <- "simulated"
output <- bind_rows(observeddata, simulateddata)
	

	# Plot the results of our model fitting.
	# biol_plot <- ggplot() +
	# 	geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) + # Observed data are points
	# 	geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
	# 	labs(x = "Time (days)", y = "Algal Biovolume")
	
	# Output the results as a ggplot2 object
	# output <- biol_plot
	return(output)
}



# SO ----------------------------------------------------------------------


# fit and plot ------------------------------------------------------------
Parameters <- c(r = 0.05, K = 10 ^ 5)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.002, K = 10 ^ 2)
UpperBound <- c(r = 5, K = 10 ^ 10) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


cells_so <- jtemp %>%
	filter(species == "SO") %>% 
	filter(total_biovolume < 1000000000) %>% 
	unite(ID, temperature, rep) %>% 
	mutate(P = "low")


all_cells_SO <-	jtemp %>%
	filter(species == "SO") %>% 
	# filter(total_biovolume < 1000000000) %>% 
	unite(Unique_ID, temperature, rep, remove = FALSE) %>% 
	rename(replicate = rep) %>% 
	select(Unique_ID, cell_density, start_time, temperature, time_since_innoc_days, replicate) %>% 
	rename(P = cell_density) %>% 
	arrange(Unique_ID) %>% 
	# filter(replicate == 1) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = Unique_ID) %>%
	arrange(days)


controldata_SO <- split(all_cells_SO, f = all_cells_SO$ID)

output_SO <- controldata_SO %>% 
	map_df(controlfit)


all_cells_SO %>% 
	filter(ID == "12_5") %>% 
	plotsinglefit(.)

jtemp %>% 
	filter(species == "SO") %>% 
ggplot(aes(x = time_since_innoc_days, y = cell_density, color = factor(temperature))) + geom_point() +
	geom_line()

output %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 10^7) %>% 
	filter(r < 2) %>% 
	# filter(temperature < 25) %>% 
	ggplot(aes( x= temperature, y = r)) + geom_point(size = 3)


write_csv(output, "data-processed/SO_fitted_rK.csv")

SO_fitted <- read_csv("data-processed/SO_fitted_rK.csv")

SO_plot <- SO_fitted %>% 
	filter(!grepl("NA", ID)) %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 10^7) %>% 
	filter(r < 2) %>% 
	# filter(temperature < 25) %>% 
	ggplot(aes( x= temperature, y = K)) + geom_point(size = 3) + theme_bw() +
	geom_smooth() +
	ylab("carrying capacity (um3/ml)") + ggtitle("Scenedemus obliquus")
	
	



# now moving onto other species here --------------------------------------

# sea <- read_csv("data-processed/sea.csv")


all_cells_TT <-	sea %>%
	filter(species == "TT") %>% 
	# filter(total_biovolume < 1000000000) %>% 
	unite(Unique_ID, temperature, rep, remove = FALSE) %>% 
	rename(replicate = rep) %>% 
	select(Unique_ID, cell_density, start_time, temperature, time_since_innoc_days, replicate, total_biovolume) %>% 
	rename(P = total_biovolume) %>% 
	arrange(Unique_ID) %>% 
	# filter(replicate == 1) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = Unique_ID) %>% 
	filter(!grepl("NA", ID)) %>% 
	arrange(days)

max(all_cells_sea$days)
ggplot(data = all_cells_sea, aes(x = days, y = P, color = ID)) + geom_point() + geom_line()

# fit and plot ------------------------------------------------------------
Parameters <- c(r = 1, K = 10 ^ 7)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, K = 10 ^ 2)
UpperBound <- c(r = 10, K = 10 ^ 9) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


controldata_TT <- split(all_cells_TT, f = all_cells_TT$ID)
str(controldata_sea)
output_TT <- controldata_TT %>% 
	map_df(controlfit)

output_TT %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	# filter(K < 50000000) %>% 
	# filter(r < 1) %>% 
	filter(temperature < 33) %>% 
	ggplot(aes(x = temperature, y = log(K))) + geom_point() + geom_smooth(method = lm) + theme_bw()

activation_energy_K <- output_TT %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 33) %>% 
	# filter(K < 50000000) %>%
	# filter(r < 1) %>%
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	do(tidy(lm(log(K) ~ inverse_temp, data = .), conf.int = TRUE)) 
write_csv(activation_energy_K, "data-processed/K_Ea.csv")

activation_energy_r <- output_TT %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	do(tidy(lm(log(r) ~ inverse_temp, data = .), conf.int = TRUE)) 






	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), r, K) %>% 
	ggplot(aes( x= temperature, y = r_mean)) + geom_point(size = 3) + theme_bw() +
	geom_errorbar(aes(ymin = r_mean - r_std.error, ymax = r_mean + r_std.error), width = 0.1) +
	geom_smooth() +
	ylab("intrinsic growth rate (r)") + ggtitle("T. tetrahele")






## something weird is going on with the 38s!
write_csv(fitted_rk_TT, "data-processed/fitted_rk_TT.csv")

all_data <- controldata_TT %>% 
	map_df(plotsinglefit, .id = "ID")

all_data %>% 
	separate(ID, into = c("temperature", "rep"), remove = FALSE) %>% 
	mutate(temperature = as.numeric(temperature)) %>% View



# Plot the results of our model fitting.
all_data <- all_data %>% 
	mutate(type = str_replace(type, "simulated", "fit"))

ggplot() +
	geom_point(data = all_data, aes(x = time, y = P, color = type)) + 
	facet_wrap( ~ ID) + theme_bw() + ylab("biovolume (um3/ml)") + xlab("days")
ggsave("figures/time_series_fits_TT.png", width = 12, height = 10)
	
	# Observed data are points
	geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
	labs(x = "Time (days)", y = "Algal Biovolume")



output_sea %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000000) %>% 
	filter(r < 1) %>% 
	filter(temperature < 33) %>% 
	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), r, K) %>% 
	ggplot(aes( x= temperature, y = r_mean)) + geom_point(size = 3) + theme_bw() +
	geom_errorbar(aes(ymin = r_mean - r_std.error, ymax = r_mean + r_std.error), width = 0.1) +
	geom_smooth() +
	ylab("intrinsic growth rate (r)") + ggtitle("T. tetrahele")



output_sea %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000000) %>% 
	filter(r < 1) %>% 
	ggplot(aes(x = temperature, y = K)) + geom_point()


# now onto CH -------------------------------------------------------------

all_cells_sea_CH <-	sea %>%
	filter(species == "CH") %>% 
	# filter(total_biovolume < 1000000000) %>% 
	unite(Unique_ID, temperature, rep, remove = FALSE) %>% 
	rename(replicate = rep) %>% 
	select(Unique_ID, cell_density, start_time, temperature, time_since_innoc_days, replicate, total_biovolume) %>% 
	rename(P = total_biovolume) %>% 
	arrange(Unique_ID) %>% 
	# filter(replicate == 1) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = Unique_ID) %>% 
	filter(!grepl("NA", ID)) %>% 
	arrange(days)

# fit and plot ------------------------------------------------------------
Parameters <- c(r = 1, K = 10 ^ 7)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, K = 10 ^ 2)
UpperBound <- c(r = 10, K = 10 ^ 9) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


controldata_sea_CH <- split(all_cells_sea_CH, f = all_cells_sea_CH$ID)
str(controldata_sea)

output_sea_CH <- controldata_sea_CH %>% 
	map_df(controlfit)

fitted_rk_CH <- output_sea_CH

all_cells_sea_CH %>% 
	filter(ID == "38_1") %>% 
	plotsinglefit(.)

CH_K_plot <- output_sea_CH %>% 
	filter(ID != "38_1") %>%  ## this one has a bad fit, taking out for now
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 500000000) %>% 
	# filter(r < 1) %>% 
	# filter(temperature < 25) %>% 
	ggplot(aes( x= temperature, y = K)) + geom_point(size = 3) + theme_bw() +
	ylab("carrying capacity, K (um3/ml)") +
	ggtitle("Chlamydomonas carrying capacity") +
	geom_smooth()
ggsave("figures/CH_K.png", width = 10, height = 8)


fitted_rk_TT %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000000) %>% 
	filter(r < 1) %>% 
	# filter(temperature < 25) %>% 
	ggplot(aes( x= temperature, y = K, color = factor(rep))) + geom_point(size = 3) + theme_bw() +
	ylab("carrying capacity, K (um3/ml)") +
	ggtitle("T. Tetrahele carrying capacity")
ggsave("figures/TT_K.png", width = 10, height = 8)

TT_K_plot <- fitted_rk_TT %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000000) %>% 
	filter(r < 1) %>% 
	filter(temperature < 35) %>% 
	ggplot(aes( x= temperature, y = K)) + geom_point(size = 3) + theme_bw() +
	ylab("carrying capacity, K (um3/ml)") +
	geom_smooth(method = lm) +
	ggtitle("T. Tetrahele carrying capacity")



all_K_combined_plot <- grid.arrange(TT_K_plot, CH_K_plot, SO_plot, nrow = 3)
ggsave("figures/all_species_combined_K.png", all_K_combined_plot, width = 8, height = 10)
?grid.arrange
