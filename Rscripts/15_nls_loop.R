sea <- read_csv("data-processed/sea_processed.csv")

TT <- sea %>% 
	filter(species == "TT") %>% 
	select(temperature, rep, cell_density, cell_volume, time_since_innoc_hours) %>% 
	mutate(time_since_innoc_hours = ifelse(is.na(time_since_innoc_hours), 12.18056, time_since_innoc_hours)) %>% 
	mutate(days = time_since_innoc_hours/24) %>% 
	unite(unique_id, temperature, rep, remove = FALSE, sep = "_")


TT_fit <- TT %>% 
	filter(cell_density != 36927) %>%
	filter(cell_density != 33992) %>% 
	filter(cell_density != 9279) %>% 
	filter(cell_density != 8924) %>% 
	filter(cell_density != 5045) %>% 
	distinct(cell_density, cell_volume, days, rep, temperature, .keep_all = TRUE)




logistic <- function(days, r, K){
	res <- K/(1 + (K/2200 - 1)*exp(-r*days))
	res
}


dat.full <- TT_fit %>% 
	rename(curve.id = unique_id)

# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
r.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
k.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits


# Loop through all curve.id.list values to estimate parameters for all curves

for(i in 1:length(curve.id.list)){
	print(i)
	
	dat<-subset(dat.full,dat.full$curve.id==curve.id.list[1])
	rvals<-seq(0,2,0.1)		
	kvals<-seq(100,10000,1000)
	mod.list<-list()
	AIC.list<-c()
	
	for(ir in 1:length(rvals)){
			for(ik in 1:length(kvals)){
			r.guess<-rvals[ir]
			k.guess<-kvals[ik]
			fit<-nlsLM(cell_density ~ k/(1 + (k/2200 - 1)*exp(-r*days)),
													 data= dat,  start=list(k = k.guess, r = r.guess),
													 lower = c(k = 100, r = 0),
													 upper = c(k = 100000, r = 2),
													 control = nls.control(maxiter=1000, minFactor=1/204800000))
				mod.list<-append(mod.list,fit)
				AIC.list<-append(AIC.list,AIC(fit))
			}
		}
	}
	
	# Identify the best model from the list and save coefficients and R^2 values
	if(!is.null(AIC.list)){
		bestmodind<-which(AIC.list==min(AIC.list))
		if(length(bestmodind)>1){
			bestmodind<-sample(bestmodind,1)
		}
		
	
		bestmod <- mod.list[[bestmodind]]
		cfs<-coef(bestmod)
		expected<-logistic(dat$days,cfs[[2]],cfs[[1]])
		rsqr<-1-sum((dat$cell_density-expected)^2)/sum((dat$cell_density-mean(dat$cell_density))^2)
	}
	
	?nlsLM
	#Save .png plot with fitted curve. File is saved with the curve.id.list value as the name
	
	#stash results		
	rsqr.list[i]<-rsqr
	K.list[i]<-cfs[[1]]
	r.list[i]<-cfs[[2]]
	n.list[i]<-length(dat$days)
}

fits2<-data.frame(curve.id.list, r.list, k.list, rsqr.list,n.list) ## constant
