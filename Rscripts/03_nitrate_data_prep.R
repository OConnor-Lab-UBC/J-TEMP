nitratejtemp <- read.csv("data-raw/nitrate_data_160104/nitratedata_160104.csv")
names(nitratejtemp)

###fit standards to absorbance curve, use curve to estimate nitrate for all samples using 'predict'
par(mfrow=c(1,1))
standards<-subset(nitratejtemp, nitratejtemp$species=="ST")
with(standards, plot(nitrate~abs))
mod<-lm(nitrate~abs, data=standards)

Ndata<-subset(nitratejtemp, nitratejtemp$species %in% c("TT", "CH"))
Ndata$nitrate<-predict(mod, data.frame(abs=Ndata$abs))

#write out porcessed data from first day
write.csv(Ndata, "data-processed/nitrate_processed.csv")
####

##repeat for redoes
nitratejtempredoes <- read.csv("data-raw/nitrate_data_160104/nitrateTTredoes_160105.csv")

#plot yesterdays standards with todays
with(standards, plot(nitrate~abs))
with(standardsredoes, points(nitrate~abs, col=2))

#rbind yesterdays standards with todays
allstandards<-rbind(standards, standardsredoes)
standardsredoes<-subset(nitratejtempredoes, nitratejtempredoes$species=="ST")
with(allstandards, plot(nitrate~abs))
mod<-lm(nitrate~abs, data=allstandards)
abline(mod)

Ndataredoes<-subset(nitratejtempredoes, nitratejtempredoes$species %in% c("TT", "CH"))
Ndataredoes$nitrate<-predict(mod, data.frame(abs=Ndataredoes$abs))

with(subset(Ndata, Ndata$species=="TT"), plot(nitrate~temp))
with(subset(Ndataredoes, Ndataredoes$species=="TT"), points(nitrate~temp, col=2))

#data from redoes match original data, so do not write over processed nitrate data
