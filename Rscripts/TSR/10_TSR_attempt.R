## jb messing around to try and figure out how adding the TSR will change the slope of the K-temp relationship

#max growth rate function
umax = function(umax_0, E_umax, k, t) umax_0 * exp(-E_umax / (k * t))

#paramter values
E_umaxn = -0.32 #activation energy 
kn = 8.62 * 10^(-5) #Boltzmann's constant

#find rate constant umax_0
tref = 273 + 15 #reference temperature
umaxref = 1000 #value of umax at reference temperature
fun = function(x) umax(x, E_umaxn, kn, tref) - umaxref
umax_0n = uniroot(fun, interval=c(0,10^12))$root #solve for umax_0

#plot umax as function of temperature
curve( umax(umax_0n, E_umaxn, kn, x), from=273.15+5, to=273.15+30, log="y", xlab="Temperature (Kelvins)", ylab="log(umax)")


MR <- function(M0R, bR, t) M0R + bR*t 
 
M0R <- 50 # base size of resource, at 5C
EP = 0.32
bR <- -0.00034882 # arbitrarily chose a negative number


## k0 is the ref K
## m is body mass
## ss is the temperature size slope (i.e. change in size due to change in temp)
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis

KMT <- function(k0, m, ss, EP, k, t) k0*((m + ss*t)^-3/4)*exp(EP/(k*t))
KMT <- function(k0, m, ss, EP, k, t) k0*(M^-3/4(exp(3em-4ep)/4kT)^-3/4)*exp(EP/(k*t))

# M^-3/4*(exp(3em-4ep)/4kT
# 				
# K <- (M0*exp(EM/kt))-3/4 * exp(EP/kt)
				
#so where em > 4/3 ep, we expect a negative k-temp relationship				
				

## draw K curve with TSR (in black)
curve(KMT(k0 = 10000, m = 100, ss = -0.10, EP = 0.32, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

## now add curve for K without TSR, shown in blue
curve(KMT(k0 = 10000, m = 100, ss = 0, EP = 0.32, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", add = TRUE, col = "blue")



#plot umax as function of temperature
curve( umax(umax_0n, E_umaxn, kn, x), from=273.15+5, to=273.15+30, xlab="Temperature (Kelvins)", ylab="log(umax)", log = "y")

 #max growth rate function
umax = function(umax_0, E_umax, k, t) umax_0 * exp(E_umax / (k * t))


?plot

#paramter values
E_umaxn = 0.32 #activation energy 
kn = 8.62 * 10^(-5) #Boltzmann's constant

#find rate constant umax_0
tref = 273 + 15 #reference temperature
umaxref = 100000 #value of umax at reference temperature
fun = function(x) umax(x, E_umaxn, kn, tref) - umaxref
umax_0n = uniroot(fun, interval=c(0,10^12))$root #solve for umax_0

#plot umax as function of temperature
curve( umax(umax_0n, E_umaxn, kn, x), from=273.15+5, to=273.15+30, xlab="Temperature (Kelvins)", ylab="log(umax)", log = "y")

curve(KMT(k0 = 1000000, M0R, br, EP, kn, x), from=273.15+5, to=273.15+30)








