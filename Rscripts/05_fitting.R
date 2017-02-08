
library(simecol)

# set up models -----------------------------------------------------------


Parameters <- c(r = 0.05, K = 10 ^ 5)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.015, K = 10 ^ 5)
UpperBound <- c(r = 2, K = 10 ^ 13) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

CRmodel <- new("odeModel",
							 main = function (time, init, parms) {
							 	with(as.list(c(init, parms)), {
							 		dp <-  r * P * (1 - (P / K))
							 		list(c(dp))
							 	})
							 },
							 parms = Parameters, # Trying Joey's empirically estimated parameters here.
							 times = c(from = 0, to = 43, by = 0.1), # the time interval over which the model will be simulated.
							 init = c(P = 2200),
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
	times(plotfittedCRmodel) <- c(from=0, to=42, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	
	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	simulateddata <- ysim
	
	# Plot the results of our model fitting.
	biol_plot <- ggplot() +
		geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) + # Observed data are points
		geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
		labs(x = "Time (days)", y = "Algal Biovolume")
	
	# Output the results as a ggplot2 object
	output <- biol_plot
	return(output)
}



# fit and plot ------------------------------------------------------------
Parameters <- c(r = 0.005, K = 10 ^ 5)

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


all_cells <-	jtemp %>%
	filter(species == "CR") %>% 
	# filter(total_biovolume < 1000000000) %>% 
	unite(Unique_ID, temperature, rep, remove = FALSE) %>% 
	mutate(P = "low") %>% 
	rename(replicate = rep) %>% 
	select(Unique_ID, cell_density, start_time, temperature, time_since_innoc_days, replicate) %>% 
	rename(P = cell_density) %>% 
	arrange(Unique_ID) %>% 
	# filter(replicate == 1) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = Unique_ID) 


controldata <- split(all_cells, f = all_cells$ID)

output <- controldata %>% 
	map_df(controlfit)


all_cells %>% 
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



# now moving onto other species here --------------------------------------

sea <- read_csv("data-processed/sea.csv")


all_cells_sea <-	sea %>%
	filter(species == "TT") %>% 
	# filter(total_biovolume < 1000000000) %>% 
	unite(Unique_ID, temperature, rep, remove = FALSE) %>% 
	rename(replicate = rep) %>% 
	select(Unique_ID, cell_density, start_time, temperature, time_since_innoc_days, replicate) %>% 
	rename(P = cell_density) %>% 
	arrange(Unique_ID) %>% 
	# filter(replicate == 1) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = Unique_ID) 


ggplot(data = all_cells_sea, aes(x = days, y = P, color = ID)) + geom_point() + geom_line()

# fit and plot ------------------------------------------------------------
Parameters <- c(r = 0.5, K = 10 ^ 5)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.2, K = 10 ^ 2)
UpperBound <- c(r = 10, K = 10 ^ 5) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


controldata_sea <- split(all_cells_sea, f = all_cells_sea$ID)

output_sea <- controldata_sea %>% 
	map_df(controlfit)


all_cells_sea %>% 
	filter(ID == "16_2") %>% 
	plotsinglefit(.)

output_sea %>% 
	separate(ID, into = c("temperature", "rep")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 50000) %>% 
	# filter(r < 2) %>% 
	# filter(temperature < 25) %>% 
	ggplot(aes( x= temperature, y = K)) + geom_point(size = 3) + scale_y_log10()

