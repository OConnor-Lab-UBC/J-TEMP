nitratejtemp <- read.csv("data-raw/nitrate_data_160104/nitratedata_160104.csv")
names(nitratejtemp)

###fit standards to absorbance curve, use curve to estimate nitrate for all samples using 'predict'
par(mfrow=c(1,1))
standards<-subset(nitratejtemp, nitratejtemp$species=="ST")
with(standards, plot(nitrate~abs))
mod<-lm(nitrate~abs, data=standards)

Ndata<-subset(nitratejtemp, nitratejtemp$species %in% c("TT", "CH"))
Ndata$nitrate<-predict(mod, data.frame(abs=Ndata$abs))
#write out porcessed data
write.csv(Ndata, "data-processed/nitrate_processed.csv")
####

