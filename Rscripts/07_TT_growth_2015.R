library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)
library(modelr)



TT <- read_csv("data-raw/TT_2015.csv")


ggplot(data = TT, aes(x = Hours.since.Innoc, y = Particles.per.ml, group = Temperature, color = factor(Temperature))) + geom_point() +
	facet_wrap( ~ N.Treatment, scales = "free") + geom_line()



TT %>% 
	filter(Temperature == 16, N.Treatment == 5) %>% 
	mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
	mutate(day = Hours.since.Innoc/24) %>% 
ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
	facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)

log(exp(1))

TT %>% 
	group_by(Temperature, N.Treatment) %>%
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>%
	ggplot(aes(x = Temperature, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
	geom_line(size = 2) + theme_bw() + facet_wrap( ~ N.Treatment) + ylab("population growth rate")

ggsave("figures/population_growth_rate_TT_2015.png")

TT_r <- TT %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 



write_csv(TT_r, "data-processed/fitted_r_TT_2015_nls.csv")

TT %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
	mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>% 
	mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
	mutate(N.Treatment = as.numeric(N.Treatment)) %>% 
	ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
	geom_line() + theme_bw() + facet_wrap( ~ Temperature)

TT %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>%
	ggplot(aes(x = Temperature, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
	geom_line(size = 2) + theme_bw() + facet_wrap( ~ N.Treatment) + ylab("population growth rate")




# now fit an activation energy --------------------------------------------

current_dataset <- TT_r %>% 
	ungroup() %>% 
	filter(N.Treatment == 2) %>% 
	rename(OriginalTraitValue = estimate) %>% 
	mutate(K = Temperature + 273.15) %>% 
	select(OriginalTraitValue, K)

current_dataset$OriginalTraitValue[current_dataset$OriginalTraitValue == 0] <- 1

## If there are negative values, substract the minimum value
MinVal <- NA
if (min(current_dataset$OriginalTraitValue)<=0){
	MinVal <- min(current_dataset$OriginalTraitValue)
	current_dataset$OriginalTraitValue <-current_dataset$OriginalTraitValue - MinVal
	current_dataset <-current_dataset[-which(current_dataset$OriginalTraitValue==0),]}



#### assign Tref as GlobalEnv
# T_ref is the standardization temperature (in K). 
# This needs to be any value below the peak of the curve.
assign("Tref", 285.15, envir = .GlobalEnv) 



# Estimate STARTING VALUES for the nls ------------------------------------

GetE <- function(tmp, rate, T.p, k=8.62e-5)
{
	# Estimate starting value for E, taking linear regression using the rise part
	# of the curve only.
	# ~~~ Parameters ~~~
	# tmp  : temperature data (in K).
	# rate : rate data corresponding to temperature above.
	# T.p  : temperature at which rate peaks, used as a cutoff point.
	# k    : Boltzmann constant.
	
	tmp.w <- which(tmp <= T.p)
	if (length(tmp.w) > 1)
	{
		m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
		return(abs(summary(m)$coefficients[2, 1]))
	} else
	{
		return(0.6)
	}
}

GetB0 <- function(tmp, rate)
{
	# Estimate starting value for the normalising constant.
	# ~~~ Parameters ~~~
	# tmp   : temperature data (in K).
	# rate  : rate data corresponding to temperature above.
	# T.ref : estimate normalising constant at this temperature (in K).
	
	if (min(tmp,na.rm=TRUE) > Tref)
	{
		return(log(min(rate[1],na.rm=TRUE)))
	} else
	{
		return(log(max(rate[which(tmp <= Tref)],na.rm=TRUE)))
	}
}


GetTpk <- function(tmp, rate)
{
	# Temperature at which the rate is maximised (estimate of T.peak).
	# ~~~ Parameters ~~~
	# tmp  : Temperature data (in K).
	# rate : Rate data corresponding to temperature above.
	
	return(max(tmp[which.max(rate)]))
}




# Schoolfield fitting -----------------------------------------------------

Schoolfield <- function(B0, E, E_D, T_h, temp) {
	
	# Boltzmann's constant. Units imply that E and E_D are in eV.
	k <- 8.62e-5
	
	# B0 is the normalization constant.    
	# E is the activation energy.
	# E_D is the de-activation energy.    
	# T_h is the temperature at which the rate-limiting enzyme 
	# is 50% active and 50% denatured due to high temperature.
	
	#     return(B0 - E/k * (1/temp - 1/Tref) - log(1 + exp((E_D/k) * (1/T_h - 1/temp)))) #Schoolfied model in original form (not with T_pk as an explicit parameter)
	
	return(B0 + log(exp((-E/k) * ((1/temp) - (1/Tref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_h - 1/temp))))) ## T_pk as an explicit parameter. FITS BETTER
	
}


B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_h_sch <- c()
T_pk_sch <- c()
P_pk_sch <- c()
AIC_sch <- c()
r_sq_sch <- c()

T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)



# schoolfield -------------------------------------------------------------

schoolfield_nls <- NA

fit_coefs <- function(current_dataset) {
schoolfield_nls <- nlsLM(
	log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
	start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
	lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
	upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
	data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024))


if(!is.na(schoolfield_nls[1])) 
{ 

	E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
	E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
	T_h_sch <- c(T_h_sch, coef(schoolfield_nls)["T_h"])
	AIC_sch<- c(AIC_sch, AIC(schoolfield_nls))
	
	# Calculate the R squared value as: 1 - (rss/tss)
	rss <- sum((exp(predict(schoolfield_nls)) - 
								current_dataset$OriginalTraitValue)^2, 
						 na.rm = TRUE)
	tss <- sum(
		(current_dataset$OriginalTraitValue - 
		 	mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
		na.rm = TRUE)
	
	if ( tss != 0 )
	{
		r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
	} else
	{
		r_sq_sch <- c(r_sq_sch, 1)
	}
	
	# Calculate the peak of the curve and its 
	# corresponding temperature value.
	curr_prediction <- predict(schoolfield_nls)
	for (j in 1:length(curr_prediction))
	{
		# If we found the maximum performance, exit the loop.
		if (curr_prediction[j] == max(curr_prediction))
		{
			break
		}
	}
	
	T_pk_sch <- c(T_pk_sch, current_dataset$K[j])

}
tidy_coefs <- tidy(schoolfield_nls)
return(tidy_coefs)
}

# Generate predictions from the model fit...
tmp_temps <- seq(min(
	floor(current_dataset$K)), 
	ceiling(max(current_dataset$K)
	), length = 200)

tmp_model <- exp(Schoolfield(
	coef(schoolfield_nls)["B0"],
	coef(schoolfield_nls)["E"],
	coef(schoolfield_nls)["E_D"],
	coef(schoolfield_nls)["T_h"],
	tmp_temps
))

ModelToPlotS <- data.frame(
	Temperature = tmp_temps - 273.15, 
	TraitValue = tmp_model
)

# Prepare the data points of the original values.
DataToPlot <- data.frame(
	Temperature = current_dataset$K - 273.15, 
	TraitValue = current_dataset$OriginalTraitValue
)
DataToPlot <- na.omit(DataToPlot)


	
ggplot(data = all_estimates, aes(x = n, y = estimate)) + geom_point(size = 3) +
	facet_wrap( ~ term, scales = "free") + theme_bw()

# plot it! ----------------------------------------------------------------

p_1 <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																											 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
																alpha = 0.7, pch = 21) + 
	geom_line(data = ModelToPlotS, 
						aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
						lwd = 1.3) +                           
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("population growth rate") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16)) +
	ggtitle("n treatment 1")

+
	ylim(0,0.05)
annotate("text", size = 3, label=             
				 	paste("R^2","sch=", sprintf("%.2f", r_sq_sch),"\nE sch=", format(coef(schoolfield_nls)["E"], digits = 3),"\nAIC sch=",format(AIC(schoolfield_nls),digits=3)), 
				 x = min(DataToPlot[, "Temperature"]),
				 y = mean(DataToPlot[, "TraitValue"]),
				 hjust=0,
				 fontface = 3)


schoolfield_nls

current_dataset <- TT_r %>% 
	ungroup() %>% 
	filter(N.Treatment != 0) %>% 
	rename(OriginalTraitValue = estimate) %>% 
	mutate(K = Temperature + 273.15) %>% 
	select(OriginalTraitValue, K, N.Treatment)

current_dataset$OriginalTraitValue[current_dataset$OriginalTraitValue <= 0] <- 0.000000000000000000001

TT_split <- split(current_dataset, f = current_dataset$N.Treatment) ## split the TT dataframe into little mini dataframes, one for each replicate - temperature combination. We will then use a map function to apply our fitting function to each dataframe individually.

coefs_TT <- TT_split %>% 
	map_df(fit_coefs, .id = "ntreatment")

write_csv(coefs_TT, "data-processed/fitted_r_TT_2015.csv")

ggplot(data = coefs_TT, aes(x = ntreatment, y = estimate)) + geom_point(size = 3) +
	facet_wrap( ~ term, scales = "free") + theme_bw()


schoolfield_nls <- NA

current_dataset <- TT_r %>% 
	ungroup() %>% 
	filter(N.Treatment == 8) %>% 
	rename(OriginalTraitValue = estimate) %>% 
	mutate(K = Temperature + 273.15) %>% 
	select(OriginalTraitValue, K, N.Treatment)

current_dataset$OriginalTraitValue[current_dataset$OriginalTraitValue <= 0] <- 0.000000000000000000001



schoolfield_nls <- nlsLM(
	log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
	start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
	lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
	upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
	data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024))


if(!is.na(schoolfield_nls[1])) 
{ 

	E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
	E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
	T_h_sch <- c(T_h_sch, coef(schoolfield_nls)["T_h"])
	AIC_sch<- c(AIC_sch, AIC(schoolfield_nls))
	
	# Calculate the R squared value as: 1 - (rss/tss)
	rss <- sum((exp(predict(schoolfield_nls)) - 
								current_dataset$OriginalTraitValue)^2, 
						 na.rm = TRUE)
	tss <- sum(
		(current_dataset$OriginalTraitValue - 
		 	mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
		na.rm = TRUE)
	
	if ( tss != 0 )
	{
		r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
	} else
	{
		r_sq_sch <- c(r_sq_sch, 1)
	}
	
	# Calculate the peak of the curve and its 
	# corresponding temperature value.
	curr_prediction <- predict(schoolfield_nls)
	for (j in 1:length(curr_prediction))
	{
		# If we found the maximum performance, exit the loop.
		if (curr_prediction[j] == max(curr_prediction))
		{
			break
		}
	}
	
	T_pk_sch <- c(T_pk_sch, current_dataset$K[j])
}
	



results <- TT %>% 
	group_by(Temperature, N.Treatment) %>%
	(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000))) 
