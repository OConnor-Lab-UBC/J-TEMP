#predict outcomes of competition between TT and CH for different temperatures and dilution rates
library(rootSolve)
library(deSolve)
library(ggplot2)
library(reshape2)

###########################################################################
#competition

#differential equations
dNNRdt = function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dN1 = N1 * ((u1 * R / (K1 + R)) - d) #consumer1
			dN2 = N2 * ((u2 * R / (K2 + R)) - d) #consumer2
			dR = d * (S - R) - Q1 * ( dN1 + d * N1) - Q2 * ( dN2 + d * N2) #resource
			return(list(c(dN1, dN2, dR)))			
		}
	)
}
	
	
#times
time = seq(from = 0, to = 100, by = 1) 

Slimiting<-27.8

umaxAF<-0.52 # at 8°C
umaxSU<-0.16 # at 8°C

KSAF<-1.6 # at 8°C
KSSU<-4.9 # at 8°C

QAF<-4.6E-7 # at 8°C
QSU<-4.6E-5 # at 8°C


state = c(N1 = 1E3, N2 = 1E3, R = Slimiting)

par(mfrow=c(1,2))

#predict for 8°C
parameters8 = c(u1 = umaxAF, K1 = KSAF, Q1 = QAF, u2 = umaxSU, K2 = KSSU, Q2 = QSU, d = 0.1, S = Slimiting)
#input N concentration in medium and original number of cells 
out = ode(y = state, times = time, func = dNNRdt, parms = parameters8)

out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column
with(subset(out.m, out.m$variable!="R"), plot(log(value, 10)  ~time, col=as.numeric(out.m$variable), xlim=c(0, 60), main="N"))#plots in base R
with(subset(out.m, out.m$variable=="R"), plot(value  ~time, col=as.numeric(out.m$variable), xlim=c(0, 60), main="R"))#plots in base R
