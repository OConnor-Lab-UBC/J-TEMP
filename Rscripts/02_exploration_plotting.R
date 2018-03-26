#### plotting the J-TEMP data

library(tidyverse)
library(plotrix)
library(broom)
# library(drc)
library(stringr)
library(lubridate)
library(growthcurve)
library(gridExtra)

jtemp <- read_csv("data-processed/Jtemp_all.csv")
nitrate_so <- read_csv("data-processed/nitrate_processed_SO.csv")
phosphate_so <- read_csv("data-processed/SO_phosphate_concentrations.csv")



jtemp %>%
	filter(species == "SO") %>% 
  filter(total_biovolume < 1000000000) %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_days, group = rep, y = total_biovolume, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature, scales = "free") + ggtitle("SO") +
	theme(axis.text.x = element_text(angle = 75, hjust = 1))


nitrate_so_2 <- nitrate_so %>% 
	mutate(new_temperature = ifelse(grepl("32", temperature), 25, temperature)) %>% 
	mutate(new_temperature = ifelse(grepl("25", temperature), 32, new_temperature)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(new_temperature = as.numeric(new_temperature)) %>% 
	rename(`nitrate concentration (uM)` = nitrate_concentration)


nitrate_plot <- ggplot(data = nitrate_so_2, aes(x = new_temperature, y = `nitrate concentration (uM)`, color = factor(rep))) + geom_point(size = 3) +
	geom_hline(yintercept = 562.42)

so_biovol_plot <- jtemp %>%
	filter(species == "SO") %>% 
	filter(start_time > "2017-01-06") %>%
	ggplot(aes(x = temperature, y = total_biovolume, color = factor(rep)), data = .) + geom_point(size = 3)




phosphate_plot <- phosphate_so %>% 
	filter(absorbance != "0.0705") %>% 
	# mutate(temperature = ifelse(species == "COMBO", 20, temperature)) %>%
	ggplot(aes(x = temperature, y = phosphate_concentration, color = factor(rep))) + geom_point(size = 3)
	



grid.arrange(so_biovol_plot, phosphate_plot, nrow = 2)




# average plots -----------------------------------------------------------

nitrate_plot <- nitrate_so_2 %>% 
	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), `nitrate concentration (uM)`) %>% 
	ggplot(data = ., aes(x = temperature, y = mean)) + geom_point(size = 3) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean+std.error), width = 0.1) +
 theme_bw() + ylab("nitrate concentration (uM)")

so_biovol_plot <- jtemp %>%
	filter(species == "SO") %>% 
	filter(start_time > "2017-01-06") %>%
	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), total_biovolume) %>% 
	ggplot(data = ., aes(x = temperature, y = mean)) + geom_point(size = 3) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean+std.error), width = 0.1) +
	theme_bw() + ylab("biovolume (um3)")



phosphate_plot <- phosphate_so %>% 
	filter(absorbance != "0.0705") %>% 
	# mutate(temperature = ifelse(species == "COMBO", 20, temperature)) %>%
	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), phosphate_concentration) %>% 
	ggplot(data = ., aes(x = temperature, y = mean)) + geom_point(size = 3) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean+std.error), width = 0.1) +
	theme_bw() + ylab("phosphate concentration (uM)")




grid.arrange(so_biovol_plot, phosphate_plot, nrow = 2)



jtemp %>% 
	filter(species == "CR", temperature != "35") %>%
	# filter(temperature < 35) %>% 
	filter(total_biovolume < 10^8) %>% 
	group_by(rep) %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, group = factor(rep), y = total_biovolume)) +
	# geom_point(size = 4) +
	geom_point() + 
	geom_line() +
	facet_wrap( ~ temperature)



# bring in the seawater species -------------------------------------------
library(tidyverse)
library(stringr)
sea_raw <- read_csv("data-processed/YangelJtemp_all.csv")
sea_raw <- read_csv("data-processed/k-temp-processed.csv")
innoc_densities <- read_csv("data-processed/jtemp_innoc_densities.csv")

sea_raw2 <- bind_rows(sea_raw, innoc_densities)

sea <- sea_raw2 %>% 
	filter(temperature != "18") %>% 
	mutate(temperature = str_replace(temperature, "24", "25")) %>% 
	mutate(temperature = as.numeric(temperature))


# plot the marine species -------------------------------------------------


sea$cell_density[sea$time_since_innoc_days < 1] <- 2200 ## just setting all the innoculation densities to the same value
write_csv(sea, "data-processed/sea2.csv")

sea %>%
	filter(temperature != "18") %>% 
	filter(species == "CH") %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_days, group = rep, y = cell_density, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature)

write_csv(sea, "data-processed/sea.csv")

# join all the data together ----------------------------------------------


all_species <- bind_rows(jtemp, sea)


# plot it all! ------------------------------------------------------------

sea %>%
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
ggplot(aes(x = time_since_innoc_days, y = cell_density, group = rep, color = factor(temperature))) +
	geom_line() + 
	geom_point() +
	facet_wrap( ~ temperature)


sea %>%
	filter(cell_density != 25415) %>% 
	filter(cell_density != 36425) %>% 
	filter(temperature != "35") %>% 
	# filter(species == "CH") %>%
	filter(time_since_innoc_days > 15) %>%
	group_by(species, temperature, rep) %>% 
	summarise(mean_density = mean(cell_density)) %>% 
	ggplot(aes(x = temperature, y = mean_density)) + geom_point() +
	geom_smooth() + facet_wrap( ~ species)

	



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

mod <- sea %>% 
	filter(species == "TT") %>%
	# filter(temperature == "5") %>% 
	group_by(temperature) %>%  
	do(tidy(x = nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .), conf.int = TRUE, conf.level = 0.95)) %>% 
	group_by(temperature) %>% 
	filter(term == "K") %>% 
	# summarise_each(funs(mean, std.error), estimate) %>% View
	ggplot(data = ., aes(x = temperature, y = estimate)) + geom_point(size = 4) +
	scale_y_log10()

## FOR SO 
mod <- jtemp %>% 
	filter(species == "SO") %>% 
	filter(total_biovolume < 1000000000) %>% 
	filter(temperature == 32) %>% 
	# group_by(rep) %>%  
	do(tidy(x = nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .), conf.int = TRUE, conf.level = 0.95)) %>% 
	group_by(temperature) %>% 
	filter(term == "K") %>% 
	# summarise_each(funs(mean, std.error), estimate) %>% View
	ggplot(data = ., aes(x = temperature, y = estimate)) + geom_point(size = 4) +
	scale_y_log10()

## 16 rep 4 can't be fit
jtemp %>% 
	filter(species == "SO") %>% 
	filter(total_biovolume < 1000000000) %>% 
	filter(temperature > 20) %>% 
	group_by(temperature, rep) %>% 
	do(tidy(nls(total_biovolume ~ Vm * time_since_innoc_days/(K+time_since_innoc_days), data = ., 
							start = list(K = max(.$time_since_innoc_days)/2, Vm = max(.$total_biovolume))))) %>% 
	filter(term == "Vm") %>%
	# group_by(temperature) %>% 
	# summarise_each(funs(mean, std.error), estimate) %>%
	# ggplot(aes(x = temperature, y = mean)) + geom_point()
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View



sea1 <- sea %>% 
	filter(time_since_innoc_days < 40)

mod5 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "5") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>% 
	# tidy(x = .,  conf.int = TRUE) %>% 
	mutate(temperature = 5)

mod8 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "8") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>%
	mutate(temperature = 8)

mod16 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "16") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>%
	mutate(temperature = 16)

mod25 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "25") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>%
	mutate(temperature = 25)

mod32 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "32") %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>%
	mutate(temperature = 32)

mod38 <- sea %>% 
	filter(species == "TT") %>%
	filter(temperature == "38") %>% 
	filter(rep != 5) %>% 
	filter(time_since_innoc_days > 0) %>% 
	group_by(rep) %>% 
	do(tidy(nls(cell_density ~ SSlogis(time_since_innoc_days, K, xmid, r), data = .))) %>%
	mutate(temperature = 38)


mods <- bind_rows(mod5, mod8, mod16, mod25, mod38)

mods %>% 
	filter(term == "K") %>% 
	filter(estimate < 50000) %>% 
ggplot(data = ., aes(x = temperature, y = estimate)) + geom_point()


mods %>% 
	filter(term == "K") %>% 
	filter(estimate < 50000) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	lm(log(estimate) ~ inverse_temp, data = .) %>% 
	tidy(conf.int = TRUE)
	summary


mods %>% 
	filter(term == "K") %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 6, alpha = 0.5, color = "#619CFF") +
	geom_smooth(method = "lm", color = "#619CFF") +
	scale_x_reverse() + xlab("temperature (1/kT)") + ylab("ln carrying capacity, K") + 
	theme_minimal() + 
	theme(axis.text.y   = element_text(size=20),
				axis.text.x   = element_text(size=20),
				axis.title.y  = element_text(size=20),
				axis.title.x  = element_text(size=20),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.ticks = element_line(size = 1)) +
	theme(panel.border = element_blank(), axis.line = element_line(colour="black", size=1, lineend="square")) +
	ylim(8.5, 11.5)



mods %>% 
	ungroup() %>% 
	filter(term == "K") %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View



# try to fit CHs with mm function -----------------------------------------

sea %>% 
	filter(species == "CH") %>% 
	group_by(temperature, rep) %>% 
	do(tidy(nls(cell_density ~ Vm * time_since_innoc_days/(K+time_since_innoc_days), data = ., 
								 start = list(K = max(.$time_since_innoc_days)/2, Vm = max(.$cell_density))))) %>% 
	filter(term == "Vm") %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	
sea %>% 
	filter(species == "CH") %>% 
	group_by(temperature, rep) %>% 
	do(tidy(nls(cell_density ~ Vm * time_since_innoc_days/(K+time_since_innoc_days), data = ., 
							start = list(K = max(.$time_since_innoc_days)/2, Vm = max(.$cell_density))))) %>% 
	filter(term == "Vm") %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 4, alpha = 0.5, color = "#619CFF") +
	geom_smooth(method = "lm", color = "#619CFF") +
	scale_x_reverse() + xlab("temperature (1/kT)") + ylab("ln carrying capacity, K") + 
	theme_minimal() + 
	theme(axis.text.y   = element_text(size=20),
				axis.text.x   = element_text(size=20),
				axis.title.y  = element_text(size=20),
				axis.title.x  = element_text(size=20),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.ticks = element_line(size = 1)) +
				theme(panel.border = element_blank(), axis.line = element_line(colour="black", size=1, lineend="square"))


?element_line

data.s <- sea %>% 
	filter(species == "CH") %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)

data.s5 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 5) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)
data.s8 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 8) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)
data.s16 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 16) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)

data.s25 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 25) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)

data.s32 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 32) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)

data.s38 <- sea %>% 
	filter(species == "CH") %>% 
	filter(temperature == 38) %>% 
	rename(density = cell_density,
				 day = time_since_innoc_days)


model.nls5 <- nls(density ~ Vm * day/(K+day), data = data.s5, 
									start = list(K = max(data.s5$day)/2, Vm = max(data.s5$density)))

model.nls8 <- nls(density ~ Vm * day/(K+day), data = data.s8, 
								 start = list(K = max(data.s8$day)/2, Vm = max(data.s8$density)))

model.nls16 <- nls(density ~ Vm * day/(K+day), data = data.s16, 
									start = list(K = max(data.s16$day)/2, Vm = max(data.s16$density)))

model.nls25 <- nls(density ~ Vm * day/(K+day), data = data.s25, 
									 start = list(K = max(data.s25$day)/2, Vm = max(data.s25$density)))

model.nls32 <- nls(density ~ Vm * day/(K+day), data = data.s32, 
									 start = list(K = max(data.s32$day)/2, Vm = max(data.s32$density)))

model.nls38 <- nls(density ~ Vm * day/(K+day), data = data.s38, 
									 start = list(K = max(data.s38$day)/2, Vm = max(data.s38$density)))


mml5 <- data.frame(day = seq(0, max(400), length.out = 50))
mml5$density <- predict(model.nls5, newdata = mml5)

mml8 <- data.frame(day = seq(0, max(400), length.out = 50))
mml8$density <- predict(model.nls8, newdata = mml8)

mml16 <- data.frame(day = seq(0, max(400), length.out = 50))
mml16$density <- predict(model.nls16, newdata = mml16)

mml25 <- data.frame(day = seq(0, max(400), length.out = 50))
mml25$density <- predict(model.nls25, newdata = mml25)

mml32 <- data.frame(day = seq(0, max(400), length.out = 50))
mml32$density <- predict(model.nls32, newdata = mml32)

mml38 <- data.frame(day = seq(0, max(400), length.out = 50))
mml38$density <- predict(model.nls38, newdata = mml38)


ggplot(data.s, aes(x = day, y = density, color = factor(temperature))) +
	theme_bw() +
	xlim(0,400) + 
	xlab("day") +
	ylab("cell density") +
	geom_point(alpha = 0.5, size = 6) +
	geom_line(data = mml5, aes(x = day, y = density), colour = "#F8766D", size = 2) +
	geom_line(data = mml8, aes(x = day, y = density), colour = "#B79F00", size = 2) +
	geom_line(data = mml16, aes(x = day, y = density), colour = "#00BA38", size = 2) +
	geom_line(data = mml25, aes(x = day, y = density), colour = "#00BFC4", size = 2) +
	geom_line(data = mml32, aes(x = day, y = density), colour = "#619CFF", size = 2) +
	geom_line(data = mml38, aes(x = day, y = density), colour = "#F564E3", size = 2)

	
library(scales)
(hue_pal()(6))



# try to fit TT with mm ---------------------------------------------------

## remove 38 rep2, all of 5?

sea %>% 
	filter(species == "TT") %>%
	filter(temperature == 25) %>%
	# filter(temperature != 38) %>% 
	# ggplot(aes(x = time_since_innoc_days, y = cell_density)) + geom_point()
	group_by(temperature, rep) %>% 
	do(tidy(nls(cell_density ~ Vm * time_since_innoc_days/(K+time_since_innoc_days), data = ., 
							start = list(K = max(.$time_since_innoc_days)/2, Vm = max(.$cell_density))))) %>% 
	filter(term == "Vm") %>% View
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	ungroup() %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View




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


