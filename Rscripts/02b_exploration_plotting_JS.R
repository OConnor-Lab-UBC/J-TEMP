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
	filter(species == "SO") %>% 
	filter(total_biovolume < 1000000000) %>% 
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
	facet_wrap( ~ temperature, scales = "free")



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
	filter(species == "CH") %>% 
	group_by(temperature, rep) %>%
	ggplot(aes(x = time_since_innoc_days, group = rep, y = cell_density, color = factor(temperature))) + geom_point(size = 4) +
	geom_line() + 
	facet_wrap( ~ temperature)

write_csv(sea, "data-processed/sea.csv")

# join all the data together ----------------------------------------------


all_species <- bind_rows(jtemp, sea)


# plot it all! ------------------------------------------------------------

names(all_species)

all_species$species<-as.factor(all_species$species)
all_species$temperature<-as.factor(all_species$temperature)
all_species<-subset(all_species, all_species$temperature %in% c("5",  "8",  "16", "25", "32","38"))
all_species<-droplevels(all_species)
summary(all_species$species)
summary(all_species$temperature)

#CH
CHdata<-subset(all_species, species %in% c("ch", "CH"))
par(mfrow=c(3,6))
for(i in 1:length(levels(CHdata$temperature))){
	subdata<-subset(CHdata, CHdata$temperature==levels(CHdata$temperature)[i])
	with(subdata, plot(cell_density~time_since_innoc_days, ylim=c(0, 50000)))
	main=levels(CHdata$temperature)[i]
}


#TT
TTdata<-subset(all_species, species %in% c("tt", "TT"))
#par(mfrow=c(2,3))
for(i in 1:length(levels(TTdata$temperature))){
	subdata<-subset(TTdata, TTdata$temperature==levels(TTdata$temperature)[i])
	with(subdata, plot(cell_density~time_since_innoc_days, ylim=c(0, 50000)))
	main=levels(TTdata$temperature)[i]
}

TTdata$time_since_innoc_days
#SO
SOdata<-subset(all_species, species %in% c("SO"))
#par(mfrow=c(2,3))
for(i in 1:length(levels(SOdata$temperature))){
	subdata<-subset(SOdata, SOdata$temperature==levels(SOdata$temperature)[i])
	with(subdata, plot(cell_density~time_since_innoc_days, ylim=c(0, 500000)))
	main=levels(SOdata$temperature)[i]
}

##########################################################
#approximate K visually so JS can present preliminary results
#For CH, take mean of data from 22 days onwards as approximation of K
##########################################################
CHK<-subset(CHdata, time_since_innoc_days>22)
CHKmean<-aggregate(CHK, by=list(CHK$temperature), FUN=mean)
CHKsd<-aggregate(CHK, by=list(CHK$temperature), FUN=sd)
par(mfrow=c(1,3))
plot(CHKmean$cell_density~levels(CHK$temperature), ylim=c(0, 40000))
segments(c(5,8,16,25,32,38), 
				 CHKmean$cell_density+CHKsd$cell_density,
				 c(5,8,16,25,32,38),
				 CHKmean$cell_density-CHKsd$cell_density)


#For TT, take mean of data from 22 days onwards as approximation of K
TTK<-subset(TTdata, time_since_innoc_days>22)
TTKmean<-aggregate(TTK, by=list(TTK$temperature), FUN=mean)
TTKsd<-aggregate(TTK, by=list(TTK$temperature), FUN=sd)
plot(TTKmean$cell_density~levels(TTK$temperature), ylim=c(0, 30000))
segments(c(5,8,16,25,32,38), 
				 TTKmean$cell_density+TTKsd$cell_density,
				 c(5,8,16,25,32,38),
				 TTKmean$cell_density-TTKsd$cell_density)

#For SO, take mean of data from 40 days onwards as approximation of K
SOK<-subset(SOdata, time_since_innoc_days>40)
SOKmean<-aggregate(SOK, by=list(SOK$temperature), FUN=mean)
SOKsd<-aggregate(SOK, by=list(SOK$temperature), FUN=sd)
plot(SOKmean$cell_density~levels(SOK$temperature), ylim=c(0, 400000))
segments(c(5,8,16,25,32,38), 
				 SOKmean$cell_density+SOKsd$cell_density,
				 c(5,8,16,25,32,38),
				 SOKmean$cell_density-SOKsd$cell_density)

names(CHK)
#Same plots but with BIOVOLUME
#For CH, take mean of data from 22 days onwards as approximation of K
CHK<-subset(CHdata, time_since_innoc_days>22)
CHKmean<-aggregate(CHK, by=list(CHK$temperature), FUN=mean)
CHKsd<-aggregate(CHK, by=list(CHK$temperature), FUN=sd)
par(mfrow=c(1,3))
plot(CHKmean$total_biovolume~levels(CHK$temperature), ylim=c(0, max(CHK$total_biovolume)/2))
segments(c(5,8,16,25,32,38), 
				 CHKmean$total_biovolume+CHKsd$total_biovolume,
				 c(5,8,16,25,32,38),
				 CHKmean$total_biovolume-CHKsd$total_biovolume)


#For TT, take mean of data from 22 days onwards as approximation of K
TTK<-subset(TTdata, time_since_innoc_days>22)
TTKmean<-aggregate(TTK, by=list(TTK$temperature), FUN=mean)
TTKsd<-aggregate(TTK, by=list(TTK$temperature), FUN=sd)
plot(TTKmean$total_biovolume~levels(TTK$temperature), ylim=c(0, max(TTK$total_biovolume)))
segments(c(5,8,16,25,32,38), 
				 TTKmean$total_biovolume+TTKsd$total_biovolume,
				 c(5,8,16,25,32,38),
				 TTKmean$total_biovolume-TTKsd$total_biovolume)

#For SO, take mean of data from 40 days onwards as approximation of K
SOK<-subset(SOdata, time_since_innoc_days>40)
SOKmean<-aggregate(SOK, by=list(SOK$temperature), FUN=mean)
SOKsd<-aggregate(SOK, by=list(SOK$temperature), FUN=sd)
plot(SOKmean$total_biovolume~levels(SOK$temperature), ylim=c(0, max(SOK$total_biovolume)))
segments(c(5,8,16,25,32,38), 
				 SOKmean$total_biovolume+SOKsd$total_biovolume,
				 c(5,8,16,25,32,38),
				 SOKmean$total_biovolume-SOKsd$total_biovolume)


