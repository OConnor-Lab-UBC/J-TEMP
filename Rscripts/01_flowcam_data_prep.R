#### flowcam file prep
#### Nov 11 2016
#### Last updated by JB



# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)


#### Step 1 #### 
## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## something like: cp **/*summary.csv /Users/Joey/Desktop/run-summaries

## step 1b, read in UniqueID key (if there is one)

# Unique_ID_key <- read.csv("/Users/Joey/Documents/chlamy-ktemp/k-temp/data-raw/PK-temp-UniqueID-key.csv")

#### Step 2: create a list of file names for each of the summaries ####

cell_files <- list.files("/Users/Joey/Documents/J-TEMP/data-raw/flowcam-summaries-nov10", full.names = TRUE)  ## find out the names of all the files in data-summary, use full.names to get the relative path for each file


names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")


#### Step 3: read in all the files!

all_cells <- map_df(cell_files, read_csv, col_names = FALSE,
										.id = "file_name")

#### Step 4: pull out just the data we want, do some renaming etc.

Jtemp <- all_cells %>% 
	rename(obs_type = X1,
				 value = X2) %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
	spread(obs_type, value) %>%
	separate(`List File`, into = c("replicate", "other"), sep = "[:punct:]") %>% 
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`,
				 cell_volume = `Volume (ABD)`)

#### Step 5: deal with the date times and separate the replicate names
Jtemp$start.time <- ymd_hms("2016-11-08 16:15:43")

Jtemp$time_since_innoc <- interval(Jtemp$start.time, Jtemp$start_time)


Jtemp_all <- Jtemp %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>%
	mutate(cell_density = as.numeric(cell_density),
				 cell_volume = as.numeric(cell_volume)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	separate(replicate, into = c("temperature", "species", "rep"), sep = c(2, 4)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	select(-other)


#### Step 6: write out the correct csv file!! yay!
write_csv(Jtemp_all, "data-processed/Jtemp_all.csv")
