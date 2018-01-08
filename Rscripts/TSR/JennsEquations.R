#necessary packages
#install.packages("rootSolve")
library(rootSolve)
library(deSolve)
library(ggplot2)
library(reshape2)

###########################################################################
#max growth rate function
umax = function(umax_0, E_umax, k, t) umax_0 * exp(-E_umax / (k * t))

#paramter values
E_umaxn = 0.65 #activation energy 
kn = 8.62 * 10^(-5) #Boltzmann's constant

#find rate constant umax_0
tref = 273 + 15 #reference temperature
umaxref = 1 #value of umax at reference temperature
fun = function(x) umax(x, E_umaxn, kn, tref) - umaxref
umax_0n = uniroot(fun, interval=c(0,10^12))$root #solve for umax_0

#plot umax as function of temperature
curve( umax(umax_0n, E_umaxn, kn, x), from=273.15+5, to=273.15+30, log="y", xlab="Temperature (Kelvins)", ylab="log(umax)")

###########################################################################
#half saturation constant function
Ks = function(Ks_0, E_Ks, k, t) Ks_0 * exp(-E_Ks / (k * t)) / (0.008 * t)

#paramter values
E_Ksn = 0.65 #activation energy 
kn = 8.62 * 10^(-5) #Boltzmann's constant

#find rate constant Ks_0
tref = 273 + 15 #reference temperature
Ksref = 1 #value of Ks at reference temperature
fun = function(x) Ks(x, E_Ksn, kn, tref) - Ksref
Ks_0n = uniroot(fun, interval=c(0,10^12))$root #solve for Ks_0

#plot Ks as function of temperature
curve( Ks(Ks_0n, E_Ksn, kn, x), from=273.15+5, to=273.15+30, log="y", xlab="Temperature (Kelvins)", ylab="log(Ks)")

###########################################################################
#Rstar function
Rstar = function(x, y, d) (x - d) / (y * d) 

#plot Rstar as a function of temperature
curve( Rstar( umax(umax_0n, E_umaxn, kn, x), Ks(Ks_0n, E_Ksn, kn, x), 0.1), from = 273.15 + 5, to = 273.15 + 30, log="y", xlab="Temperature (Kelvins)", ylab="log(Rstar)")

###########################################################################
#numerically solve and plot one species temporal dynamics

#times
time = seq(from = 0, to = 100, by = 0.01) 

#parameters
temp = 273.15 + 5 #temperature
parameters = c(u = umax(umax_0n, E_umaxn, kn, temp), K = Ks(Ks_0n, E_Ksn, kn, temp), d = 0.1, S = 2, Q = 0.1)

#initial state
state = c(N = 0.1, R = 10)

#differential equations
dNRdt = function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dN = N * (u / (K + R) * R - d) 			#consumer
			dR = d * (S - R) - Q *( dN + d * N) 	#resource
			return(list(c(dN,dR)))			
		}
	)
}
	
out = ode(y = state, times = time, func = dNRdt, parms = parameters)
plot(out, xlab="time", ylim=c(0,40))

###########################################################################
#competition

#times
time = seq(from = 0, to = 100, by = 0.01) 

#parameters
temp = 273.15 + 5 #temperature
parameters = c(u1 = umax(umax_0n, E_umaxn, kn, temp), K1 = Ks(Ks_0n, E_Ksn, kn, temp), Q1 = 0.1, u2 = umax(umax_0n, E_umaxn+0.01, kn, temp), K2 = Ks(Ks_0n, E_Ksn, kn, temp), Q2 = 0.1,d = 0.1, S = 2)

#initial state
state = c(N1 = 4, N2 = 4, R = 10)

#differential equations
dNNRdt = function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dN1 = N1 * (u1 / (K1 + R) * R - d) #consumer1
			dN2 = N2 * (u2 / (K2 + R) * R - d) #consumer2
			dR = d * (S - R) - Q1 * ( dN1 + d * N1) - Q2 * ( dN2 + d * N2) #resource
			return(list(c(dN1, dN2, dR)))			
		}
	)
}
	
out = ode(y = state, times = time, func = dNNRdt, parms = parameters)
#plot(out, xlab="time", ylim=c(0,40))

out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column
    
p <- ggplot(out.m, aes(time, value, color = variable)) + geom_point()
print(p)


 
	