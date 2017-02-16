setwd("~/Dropbox/interactions and warming/summer 2015 data/T. tetra/Aug 24 Scott/") # change to where you saved this directory

A10data<-read.csv(file="TT week of aug 24.csv")

A10data <- TT

A10data$N.Treatment<-as.factor(A10data $N.Treatment)
A10data$Temperature <-as.factor(A10data $Temperature)


A10exp<-A10data



###############################
par(mfrow=c(2,3))
########
########
########
#13
rbyN13<-data.frame(Nconc=c(0:8), Temp=seq(13, 13, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))

with(subset(A10exp, A10exp$Temperature=="13"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", xlab="Hours since inocculation", las=1, cex.axis=0.8)) 

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("13°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp$N==levels(A10exp$N.Treatment)[i] & A10exp$Temperature==13)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata $Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN13[i,4] <-coef(curve)[2]
rbyN13[i,3] <-coef(curve)[1]
}


#16
rbyN16<-data.frame(Nconc=c(0:8), Temp=seq(16, 16, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))

with(subset(A10data, A10data$Temperature=="16"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", xlab="Hours since inocculation", las=1, cex.axis=0.8)) 

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("16°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp $N==levels(A10exp $N.Treatment)[i] & A10exp $Temperature==16)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata $Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN16[i,4] <-coef(curve)[2]
rbyN16[i,3] <-coef(curve)[1]
}


#19

rbyN19<-data.frame(Nconc=c(0:8), Temp=seq(19, 19, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))

with(subset(A10data, A10data$Temperature=="19"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", xlab="Hours since inocculation", las=1, cex.axis=0.8)) 

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("19°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp $N==levels(A10exp $N.Treatment)[i] & A10exp $Temperature==19)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata$Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN19[i,4] <-coef(curve)[2]
rbyN19[i,3] <-coef(curve)[1]
}

#22
rbyN22<-data.frame(Nconc=c(0:8), Temp=seq(22, 22, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))


with(subset(A10data, A10data$Temperature=="22"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", xlab="Hours since inocculation", las=1, cex.axis=0.8)) 

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("22°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp $N==levels(A10exp $N.Treatment)[i] & A10exp $Temperature==22)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata $Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN22[i,4] <-coef(curve)[2]
rbyN22[i,3] <-coef(curve)[1]
}

#25
rbyN25<-data.frame(Nconc=c(0:8), Temp=seq(25, 25, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))

with(subset(A10data, A10data$Temperature=="25"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", xlab="Hours since inocculation", las=1, cex.axis=0.8)) 

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("25°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp $N==levels(A10exp $N.Treatment)[i] & A10exp $Temperature==25)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata $Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN25[i,4] <-coef(curve)[2]
rbyN25[i,3] <-coef(curve)[1]
}

#28
rbyN28<-data.frame(Nconc=c(0:8), Temp=seq(28, 28, length.out=9), r=seq(0, 0, length.out=9), a=seq(0, 0, length.out=9))


with(subset(A10data, A10data$Temperature=="28"), plot(Particles.per.ml~ Hours.since.Innoc, col=as.numeric(N.Treatment), pch=as.numeric(N.Treatment), ylab="cells per ml", ylim=c(0, 30000), xlab="Hours since inocculation", las=1, cex.axis=0.8))  

legend("topleft", title="N dilution", levels(as.factor(A10exp $N.Treatment)), col=1:10, pch=1:10, lty=1, bty="n", cex=0.8)
mtext("28°C", 3, -2)

a1 = 0.01
b1 = 100

for(i in c(1:9)){
		
	curvedata<-subset(A10exp, A10exp $N==levels(A10exp $N.Treatment)[i] & A10exp $Temperature==28)
	
	curve = nls(Particles.per.ml~ 75 * (1+a)^(Hours.since.Innoc), data= curvedata,  start=list(a=a1), control = nls.control(maxiter=100, minFactor=1/204800000))

a1 = coef(curve)[1]
b1 = coef(curve)[2]

	#plot the line 
new = data.frame(Hours.since.Innoc = seq(min(curvedata$Hours.since.Innoc),max(curvedata $Hours.since.Innoc),len=200))#make dummy x data
lines(new$Hours.since.Innoc,predict(curve, newdata=new), col=i, lty=1)#draw predicted line connecting each point

rbyN28[i,4] <-coef(curve)[2]
rbyN28[i,3] <-coef(curve)[1]
}





