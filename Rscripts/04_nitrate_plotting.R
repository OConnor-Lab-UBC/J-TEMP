Ndata <- read.csv("data-processed/nitrate_processed.csv")

par(mfrow=c(1,2))
TTNdata<-subset(Ndata, Ndata$species=="TT")
CHNdata<-subset(Ndata, Ndata$species=="CH")
with(TTNdata, plot(abs~temp))
with(CHNdata, plot(abs~temp))


pdf("R-star_across_T.pdf", 6, 5)
par(mfrow=c(2,2))
par(mar=c(3,4,1, 0.5), oma=c(1, 0, 1, 0))
with(TTNdata, plot(nitrate~temp, ylim=c(0,15), xlab="", ylab="Nitrate remaining, umol", main="TT"))
with(CHNdata, plot(nitrate~temp, ylim=c(0,15), xlab="", ylab="",  main="CH"))

plot(CHKmean$cell_density~levels(CHK$temperature), ylim=c(0, 40000), xlab="Temperature", ylab="Carrying capacity, cells/ml")
segments(c(5,8,16,25,32,38), 
				 CHKmean$cell_density+CHKsd$cell_density,
				 c(5,8,16,25,32,38),
				 CHKmean$cell_density-CHKsd$cell_density)


#For TT, take mean of data from 22 days onwards as approximation of K
plot(TTKmean$cell_density~levels(TTK$temperature), ylim=c(0, 30000), xlab="Temperature", ylab="")
segments(c(5,8,16,25,32,38), 
				 TTKmean$cell_density+TTKsd$cell_density,
				 c(5,8,16,25,32,38),
				 TTKmean$cell_density-TTKsd$cell_density)
dev.off()
