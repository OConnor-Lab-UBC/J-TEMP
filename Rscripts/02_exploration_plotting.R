#### plotting the J-TEMP data

library(tidyverse)
library(plotrix)
library(broom)
# library(drc)
library(stringr)
library(lubridate)
library(growthcurve)

jtemp <- read_csv("data-processed/Jtemp_CR_all.csv")

jtemp %>%
	filter(species == "CR", temperature == 38) %>% View
	arrange()
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_hours, group = rep, y = total_biovolume, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature) + ggtitle("Chlamydomonas reinhardtii")

jtemp %>% 
	filter(species == "CR", temperature != "35") %>%
	group_by(rep) %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, group = factor(rep), y = cell_density)) +
	# geom_point(size = 4) +
	geom_point() + 
	geom_line() +
	facet_wrap( ~ temperature)



# bring in the seawater species -------------------------------------------

sea_raw <- read_csv("data-processed/YangelJtemp_all.csv")
innoc_densities <- read_csv("data-processed/jtemp_innoc_densities.csv")

sea_raw2 <- bind_rows(sea_raw, innoc_densities)



unique(sea_raw$temperature)

sea <- sea_raw2 %>% 
	filter(temperature != "18") %>% 
	mutate(temperature = str_replace(temperature, "24", "25")) %>% 
	mutate(temperature = as.numeric(temperature))


# plot the marine species -------------------------------------------------


sea$cell_density[sea$time_since_innoc_days < 1] <- 2200 ## just setting all the innoculation densities to the same value


sea %>%
	filter(temperature != "18") %>% 
	filter(species == "TT") %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_days, group = rep, y = cell_density, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature)

write_csv(sea, "data-processed/sea.csv")

# join all the data together ----------------------------------------------


all_species <- bind_rows(jtemp, sea)


### let's go in and fix the time zero cell densities

sea2 <- sea %>% 
	mutate(cell_density = ifelse(time_since_innoc_days < 1, 2200, cell_density)) %>%
	arrange(time_since_innoc_days)
	


# plot it all! ------------------------------------------------------------

sea2 %>%
	filter(cell_density != 25415) %>% 
	filter(cell_density != 36425) %>% 
	# filter(species %in% c("CH", "TT", "CR"), temperature != 35) %>% 
	filter(temperature != "35") %>% 
	filter(species == "CH") %>% 
	mutate(month = month(start_time)) %>%
	mutate(day = day(start_time)) %>% 
	mutate(year = 2016) %>% 
	unite(month_day, month, day, year) %>%
	mutate(month_day = mdy(month_day)) %>% 
	group_by(rep) %>%
ggplot(aes(x = time_since_innoc_hours, y = cell_density, group = rep, color = factor(temperature))) +
	geom_line() + 
	geom_point() +
	facet_wrap( ~ temperature)
	



all_species %>%
	filter(cell_density != 25415) %>% 
	filter(cell_density != 36425) %>% 
	# filter(species %in% c("CH", "TT", "CR"), temperature != 35) %>% 
	filter(temperature != "35") %>% 
	filter(species == "CR") %>% 
	mutate(month = month(start_time)) %>%
	mutate(day = day(start_time)) %>% 
	mutate(year = 2016) %>% 
	unite(month_day, month, day, year) %>%
	mutate(month_day = mdy(month_day)) %>% 
	group_by(temperature, month_day, species) %>% 
	summarise_each(funs(mean, std.error), total_biovolume) %>% 
	ggplot(aes(x = month_day, y = mean, group = temperature, color = factor(temperature))) +
geom_point(size = 5) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 1) +
	geom_line() + 
	facet_wrap( ~ temperature) + 
	xlab("date") + ylab("mean cell density") + 
	theme_minimal() +
	theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_rect(colour = "grey", size=1))




# estimate and plot r and K -----------------------------------------------

all_species %>% 
	filter(species == "CH") %>%
	filter(temperature == "38") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(temperature) %>%  
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>% 
	group_by(temperature) %>% 
	filter(term == "K") %>% View
	# summarise_each(funs(mean, std.error), estimate) %>% View
	ggplot(data = ., aes(x = temperature, y = estimate)) + geom_point(size = 4) +
	scale_y_log10()



+ facet_wrap( ~ term, scales = "free") +
	geom_errorbar(aes(ymin = estimate - std.error, ymax= estimate + std.error))



nls(norm_div~K/(1+exp(r*norm_time)), data=df, 
		start=list(K=1, r=-1.6), 
		control=list(maxiter=1000, minFactor=.00000000001))



# trying again! -----------------------------------------------------------

CR <- all_species %>% 
	filter(species == "CR") %>%
	filter(time_since_innoc_days > 1)

	mutate(time_since_innoc_days = ifelse(time_since_innoc_days < 0, 0, time_since_innoc_days))

CR_innoc <- read_csv("data-raw/innoc_CR.csv") ### adding innoc data

CR_all <- bind_rows(CR, CR_innoc)	

CR_dat <- CR_all %>% 
	dplyr::select(time_since_innoc_days, cell_density, total_biovolume, temperature, rep)


CR_dat %>%
	filter(cell_density != 25415) %>% 
	filter(cell_density != 36425) %>% 
	# filter(species %in% c("CH", "TT", "CR"), temperature != 35) %>% 
	filter(temperature != "35") %>% 
	group_by(rep) %>% 
	ggplot(aes(x = time_since_innoc_days, y = total_biovolume, group = rep, color = factor(temperature))) +
	geom_point(size = 5)+
	geom_line() + facet_wrap( ~ temperature, scales = "free") +
	scale_y_log10()


CR_mod <- CR_dat %>% 
	filter(temperature == "25", time_since_innoc_days < 5) %>% 
	nls(cell_density ~ 
				SSlogis(time_since_innoc_days, Asym, xmid, scal), 
			data = .)


newstart <- list(Asym = max(CR_dat$cell_density),
								 xmid = mean(CR_dat$time_since_innoc_days),
								 scal=1)
CE.fullmod <- nls(cell_density ~ SSlogis(time_since_innoc_days), data=CR_dat,
										 start=newstart)

summary(CE.fullmod)


summary(CR_mod)


### try something else

###Log fit - be sure to use quotes around the variable names in the call
log.fit <- function(dep, ind, yourdata){
	#Self-starting...
	
	y <- yourdata[, dep]
	x <- yourdata[, ind]
	
	log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
	
	#C
	C <- summary(log.ss)$coef[1]
	#a
	A <- exp((summary(log.ss)$coef[2]) * (1/summary(log.ss)$coef[3]))
	#k
	K <- (1 / summary(log.ss)$coef[3])
	
	plot(y ~ x, main = "Logistic Function", xlab=ind, ylab=dep)
	lines(0:max(x), predict(log.ss, data.frame(x=0:max(x))), col="red")
	
	r1 <- sum((x - mean(x))^2)
	r2 <- sum(residuals(log.ss)^2)
	
	r_sq <- (r1 - r2) / r1
	
	out <- data.frame(cbind(c(C=C, a=A, k=K, R.value=sqrt(r_sq))))
	names(out)[1] <- "Logistic Curve"
	
	return(out)
}


cr_sub <- CR_dat %>% 
	filter(temperature == "16", time_since_innoc_days < 13) 

log.fit(dep = "cell_density", ind = "time_since_innoc_days", yourdata = cr_sub)

str(cr_sub)



# trying something else yet again! this time with growth curve ------------

all_species %>% 
	filter(species == "CH") %>%
	ggplot(aes(x = time_since_innoc_days, y = cell_density, group = rep)) + geom_point() +
	facet_wrap( ~ temperature) + geom_line()


fits8 <- all_species %>% 
	filter(species == "TT", temperature == "8") %>%
	fit_growth(time_since_innoc_days, cell_density, type = "logistic")

K5 <- as.data.frame(fits5$parameters$max_growth) %>% 
	t(.)
K16 <- as.data.frame(fits16$parameters$max_growth)
K25 <- as.data.frame(fits25$parameters$max_growth)
K32 <- as.data.frame(fits32$parameters$max_growth)
K38 <- as.data.frame(fits38$parameters$max_growth)


allK <- bind_rows(K5, K16, K25, K32, K38)


?fit_growth

str(fits)

plot(fits, show_raw=TRUE, show_maxrate=FALSE, show_asymptote=FALSE)

autoplot(fits)

str(all_species)

all_species %>% 
	filter(species == "TT", temperature == "25") %>% 
	ggplot(aes(x = time_since_innoc_days, y = cell_density)) + geom_point() +
	stat_growthcurve(type = "logistic")

?stat_growthcurve
