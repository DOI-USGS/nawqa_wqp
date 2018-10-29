# GAP analysis script for the water quality coordination effort
# this script filters the WQ data based on filtering criteria

library(tidyverse)
library(feather)

#***********************
# batch mode or manually
runMode<-'batch'
if (runMode=='batch'){
  setwd('D:/abock/NAWQA/WQR/nawqa_wqp') # Set locally for running in batch mode for quick processing of all states
  args<-commandArgs(trailingOnly = T)
  state<-args[1]
}else{
  state<-'Alabama'
}

# WQ constituents to consider in this script
Cons<-c("Temp","DO","pH","SC")

print(state)
# State abbreviation and data files
state2<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state2, state.name)][1]
data_files <- unlist(lapply(Cons,function (con) {dir(paste('4_gap_analysis/in/Cons/',con,sep=""), pattern=paste(stateAb,'_',con, sep=''),full.names = T)})) %>%
  grep('feather', ., value=TRUE)

# Filtering methods
# combine WQ data into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    cname <- unlist(strsplit(dfile,"/"))[4]
    print(cname)
    dat <- feather::read_feather(dfile) %>%
    mutate(ConstituentGroup=cname) %>%
      # Subset by records that are quality control
      filter(!grepl("Quality Control",ActivityTypeCode)) %>%
      # Fill in with ActivityStartDate in case ActivityStartDateTime has no value
      mutate(ActivityStartDateTime=if_else(is.na(ActivityStartDateTime),
                                           as.POSIXct(paste(ActivityStartDate,"12:00:00",sep=" "),tz="UTC"),
                                           ActivityStartDateTime)) %>%
    # Filter to just the dates we are interested in, this will cut down number of sites
    filter(ActivityStartDateTime > "2010-01-01" & ActivityStartDateTime < "2016-01-01") %>%
    # Aggregate to daily
    group_by(OrganizationIdentifier,MonitoringLocationIdentifier,ConstituentGroup,ActivityStartDateTime) %>%
    summarize(measure=mean(ResultMeasureValue),totalObs=n())
  }))

# count the number of target constituents per sample
constituents_per_sample <- combined_data %>%
  mutate(Year = as.integer(format(ActivityStartDateTime, '%Y'))) %>%
  # Group by location, date, and year
  group_by(MonitoringLocationIdentifier, ActivityStartDateTime, Year) %>%
  # how many Cons per sample
  summarize(NumConstituentsPerSample=length(unique(ConstituentGroup)))

print(unique(constituents_per_sample$NumConstituentsPerSample))

# Designate time stamps for files
# Year, month, week, day, and quarter
TimeAtt_perSample<- constituents_per_sample %>%
  select(MonitoringLocationIdentifier,ActivityStartDateTime,NumConstituentsPerSample) %>%
  mutate(Year = as.integer(format(ActivityStartDateTime, '%Y'))) %>%
  mutate(Month = as.integer(format(ActivityStartDateTime, '%m'))) %>%
  mutate(Week = as.integer(strftime(ActivityStartDateTime, format = "%V"))) %>%
  mutate(Day = as.integer(format(ActivityStartDateTime, '%d'))) %>%
  arrange(MonitoringLocationIdentifier,ActivityStartDateTime) %>%
  mutate(NumQtrs1=ifelse(Month %in% c(4:6),2,ifelse(Month %in% c(7:9),3,ifelse(Month %in% c(10:12),4,1)))) %>%
  mutate(NumQtrs2=ifelse(Month %in% c(5:7),2,ifelse(Month %in% c(8:10),3,ifelse(Month %in% c(1,11,12),4,1)))) %>%
  mutate(NumQtrs3=ifelse(Month %in% c(6:8),2,ifelse(Month %in% c(9:11),3,ifelse(Month %in% c(1,2,12),4,1)))) %>%
  group_by(MonitoringLocationIdentifier) %>%
  # Calculate the number of days between sampling dates
  mutate(diff=ActivityStartDateTime-lag(ActivityStartDateTime),
         diff_days=as.numeric(diff,units="days"),
         diff_weeks=ifelse(Week<5,(Week+52)-lag(Week),Week-lag(Week))) %>%
  select(-diff) %>% replace_na(list(diff_days=0,diff_weeks=0))

# Write out summaries to the WQ folder
write_feather(TimeAtt_perSample,paste('4_gap_analysis/in/Cons/WQ/',stateAb[1],'_samples.feather',sep=''))
write.csv(TimeAtt_perSample,paste('4_gap_analysis/in/Cons/WQ/',stateAb[1],'_samples.csv',sep=''),row.names=F,quote=F)

# Summarize sampling time series for each site
siteSummary<-TimeAtt_perSample %>% group_by(MonitoringLocationIdentifier,NumConstituentsPerSample,Year) %>%
  summarize(numQtrs=max(c(NumQtrs1,NumQtrs2,NumQtrs3)),numMonths=n_distinct(Month),numWeeks=n_distinct(Week),
            numWeekSamp=sum(diff_weeks ==0 | diff_weeks < 6))%>% ungroup() %>%
  group_by(MonitoringLocationIdentifier)%>%
  summarize(minSamps=min(NumConstituentsPerSample),maxSamps=max(NumConstituentsPerSample),
            numQtrs=sum(numQtrs),numMonths=sum(numMonths),
            numWeeks=sum(numWeekSamp),numYrs=n_distinct(Year),lastYear=max(Year))

# Site Criteria
# number of years, number of weeks in year, number of days in year
# Level 1 - daily data (>5 years)
# Level 2 - weekly (>5 years)
# Level 3 - bi-weekly (>5 years)
# Level 4 - monthly (>5 years)
# Level 5 - seasonal (> 5 years)
# Level a - same as 1-5 above, but < 3 cons
# Level b - same as 1-5 above, but < 3 cons, but missing starting years
# Level c - same as 1-5 above,  but < 3 cons, and missing 75% of sampling dates
# Level d - 4 cons, 50-100% of sampling frequency criteria met
Tiered_sites<-siteSummary %>% group_by(MonitoringLocationIdentifier) %>%
  summarize(Tier=case_when(minSamps >3 & numYrs >= 5 & numWeeks >= 260 ~ "Level 2a",
                           minSamps >3 & numYrs >= 5 & numWeeks >= 130 ~ "Level 3a",
                           minSamps >3 & numYrs >= 5 & numMonths >= 48 ~ "Level 4a",
                           minSamps >3 & numYrs >= 5 & numQtrs >= 20 ~ "Level 5a",
                           # second set
                           minSamps >=3 & numYrs >= 5 & numWeeks >= 260 ~ "Level 2b",
                           minSamps >=3 & numYrs >= 5 & numWeeks >= 130 ~ "Level 3b",
                           minSamps >=3 & numYrs >= 5 & numMonths >= 60 ~ "Level 4b",
                           minSamps >=3 & numYrs >= 5 & numQtrs >= 20 ~ "Level 5b",
                           # third set
                           minSamps >3 & numYrs < 5 & numWeeks >= (numYrs*52) & max(lastYear)==2015 ~ "Level 2c",
                           minSamps >3 & numYrs < 5 & numWeeks >= (numYrs*26) & max(lastYear)==2015 ~ "Level 3c",
                           minSamps >3 & numYrs < 5 & numMonths >= (numYrs*12)  & max(lastYear)==2015 ~ "Level 4c",
                           minSamps >3 & numYrs < 5 & numQtrs >= (numYrs*4) & max(lastYear)==2015 ~ "Level 5c",
                           # fourth set
                           minSamps >3 & numYrs >= 3 & numWeeks >= 130 & numWeeks < 260 ~ "Level 2d",
                           minSamps >3 & numYrs >= 3 & numWeeks >= 65 & numWeeks < 130 ~ "Level 3d",
                           minSamps >3 & numYrs >= 3 & numMonths >= 24 & numMonths < 48 ~ "Level 4d",
                           minSamps >3 & numYrs >= 3 & numQtrs >= 10 & numQtrs < 20 ~ "Level 5d"))#,

Tiered_sites %>% group_by(Tier) %>% summarize(n=n())

# Write out results
write_feather(Tiered_sites,paste('4_gap_analysis/in/Tiers/WQ/',stateAb[1],'_Tiers.feather',sep=''))
write.csv(Tiered_sites,paste('4_gap_analysis/in/Tiers/WQ/',stateAb[1],'_Tiers.csv',sep=''),row.names=F,quote=F)

# need to polish up this next part for both Nut/WQ/SingleCon

# An example of how to join with the site data to prep for the map
site_files <- unlist(lapply(Cons,function (con) {dir(paste('4_gap_analysis/in/Site/',con,sep=""), pattern=paste(stateAb),full.names = T)})) %>%
  grep('feather', ., value=TRUE)
# Combine Sites
combined_sites <- bind_rows(lapply(site_files, function(sfile) {dat <- feather::read_feather(sfile) %>%
 mutate_at(vars(starts_with('StateCode')),funs(as.character))}))

# unique list of sites
combined_sites <- combined_sites %>% distinct(MonitoringLocationIdentifier,.keep_all=T)
combined_sites<-inner_join(x=combined_sites,y=Tiered_sites,by="MonitoringLocationIdentifier")

# Write out results
write_feather(combined_sites,paste('4_gap_analysis/in/Site/WQ/',stateAb[1],'_sites.feather',sep=''))
write.csv(combined_sites,paste('4_gap_analysis/in/Site/WQ/',stateAb[1],'.csv',sep=''),row.names=F,quote=F)




