# GAP analysis script for the water quality coordination effort

library(tidyverse)
library(feather) # install with install.packages('feather')
library(scipiper) # install with devtools::install_github('USGS-R/scipiper')
library(XLConnect) # read xlsx files

#***********************
# batch mode or manually
runMode<-'man'
if (runMode=='batch'){
  setwd('D:/abock/NAWQA/WQR/nawqa_wqp') # Set locally for running in batch mode for quick processing of all states
  args<-commandArgs(trailingOnly = T)
  state<-args[1]
}else{
  state<-'Alabama'
}

# State abbreviation and data files
state2<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state2, state.name)]
data_files <- dir('4_gap_analysis/in/Cons', pattern=paste(stateAb,'_*', sep=''),full.names = T) %>%
  grep('feather', ., value=TRUE)

# Filtering methods
# combine N/P data into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {dat <- feather::read_feather(dfile)})) %>%
  # Subset by records that are quality control
  filter(!grepl("Quality Control",ActivityTypeCode))

# Grep on gap names with dissolved/total p
TP<-combined_data %>% filter(grepl("\\bTotal\\b",GAP_name)) %>% mutate(ConstituentGroup="Total Phosphorus")
DP<-combined_data %>% filter(grepl("\\bDissolved\\b",GAP_name)) %>% mutate(ConstituentGroup="Dissolved Phosphorus")
Nit<-combined_data %>% filter(ConstituentGroup=="nitrate")

print(dim(combined_data))
combined_data<-rbind(TP,DP,Nit)

# count the number of target constituents per sample
constituents_per_sample <- combined_data %>%
  select(ConstituentGroup,MonitoringLocationIdentifier,ResultMeasureValue,ActivityStartDateTime) %>%
  mutate(Year = as.integer(format(ActivityStartDateTime, '%Y'))) %>%
  # Group by location, date, and year
  group_by(MonitoringLocationIdentifier, ActivityStartDateTime, Year) %>%
  # how many Cons per sample
  summarize(NumConstituentsPerSample=length(unique(ConstituentGroup)))

print(unique(constituents_per_sample$NumConstituentsPerSample))

# Write out summaries
write_feather(constituents_per_sample,paste('4_gap_analysis/in/Cons/',stateAb,'_samples.feather',sep=''))
write.csv(constituents_per_sample,paste('4_gap_analysis/in/Cons/',stateAb,'_samples.csv',sep=''),row.names=F,quote=F)

# Seasonal filtering
sea_samples<-constituents_per_sample %>%
  mutate(Month=as.integer(format(ActivityStartDateTime,'%m'))) %>%
  # Designate three possible seasonal sampling schemes
  mutate(Qtr=ifelse(Month %in% c(4:6),2,ifelse(Month %in% c(7:9),3,ifelse(Month %in% c(10:12),4,1)))) %>%
  mutate(Qtr2=ifelse(Month %in% c(5:7),2,ifelse(Month %in% c(8:10),3,ifelse(Month %in% c(1,11,12),4,1)))) %>%
  mutate(Qtr3=ifelse(Month %in% c(6:8),2,ifelse(Month %in% c(9:11),3,ifelse(Month %in% c(1,2,12),4,1)))) %>%
  group_by(MonitoringLocationIdentifier,Year) %>%
  mutate(NumQtrs1=n_distinct(Qtr),NumQtrs2=n_distinct(Qtr2),NumQtrs3=n_distinct(Qtr3))

# Write out results
write_feather(sea_samples,paste('4_gap_analysis/in/Tiers/',stateAb,'_samples.feather',sep=''))
write.csv(sea_samples,paste('4_gap_analysis/in/Tiers/',stateAb,'_samples.csv',sep=''),row.names=F,quote=F)




