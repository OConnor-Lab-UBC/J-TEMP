


sea <- read_csv("data-processed/sea_processed2.csv")


sea %>% 
	filter(temperature == 5) %>% 
	filter(time_since_innoc_days > 30) %>% 
	summarise(mean_size = mean(cell_volume))


sea %>% 
	filter(temperature < 32) %>% 
	filter(time_since_innoc_days > 25) %>% 
	group_by(rep, temperature) %>% 
	summarise(mean_size = mean(cell_volume)) %>% 
	# ungroup() %>% 
	# ggplot(aes(x = temperature, y = mean_size)) +geom_point() +
	# geom_smooth(method = "lm")
	lm(mean_size ~ temperature, data = .) %>% summary()

sea %>% 
	filter(temperature < 32) 



(-15.991 /1173)*100

