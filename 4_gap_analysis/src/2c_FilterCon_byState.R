# GAP analysis script for the water quality coordination effort
# this script filters the data for a single constituent based on filtering criteria

# To Do - create files for other P and N types

library(tidyverse)
library(feather)

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

#Possible cons<-c("TP","DP","N","Temp","DO","pH","SC")
Con<-"SC"

print(state)
# State abbreviation and data files
state2<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state2, state.name)][1]
data_files <- dir(paste('4_gap_analysis/in/Cons/',Con,sep=""), pattern=paste(stateAb,'_*', sep=''),full.names = T) %>%
  grep(paste0(Con,'.feather'), ., value=TRUE)

# Keeping the function for now, even though there is just one file for each CON
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {dat <- feather::read_feather(dfile)})) %>%
  # Subset by records that are quality control
  filter(!grepl("Quality Control",ActivityTypeCode)) %>%
  # Fill in with ActivityStartDate in case ActivityStartDateTime has no value
  mutate(ActivityStartDateTime=if_else(is.na(ActivityStartDateTime),
                                       as.POSIXct(paste(ActivityStartDate,"12:00:00",sep=" "),tz="UTC"),
                                       ActivityStartDateTime)) %>%
  filter(ActivityStartDateTime > "2010-01-01" & ActivityStartDateTime < "2016-01-01") %>%
  group_by(OrganizationIdentifier,MonitoringLocationIdentifier,ActivityStartDateTime) %>%
  summarize(measure=mean(ResultMeasureValue),totalObs=n())

# Designate time stamps for files
# Year, month, week, day, and quarter
TimeAtt_perSample<- combined_data %>%
  select(MonitoringLocationIdentifier,ActivityStartDateTime) %>%
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


# Write out summaries
write_feather(TimeAtt_perSample,paste('4_gap_analysis/in/Cons/',Con,'/',stateAb[1],'_samples.feather',sep=''))
write.csv(TimeAtt_perSample,paste('4_gap_analysis/in/Cons/',Con,'/',stateAb[1],'_samples.csv',sep=''),row.names=F,quote=F)

# Summarize sampling time series for each site
siteSummary<-TimeAtt_perSample %>% group_by(MonitoringLocationIdentifier,Year) %>%
  summarize(numQtrs=max(c(NumQtrs1,NumQtrs2,NumQtrs3)),numMonths=n_distinct(Month),numWeeks=n_distinct(Week),
            numWeekSamp=sum(diff_weeks ==0 | diff_weeks < 6))%>% ungroup() %>%
  group_by(MonitoringLocationIdentifier)%>%
  summarize(numQtrs=sum(numQtrs),numMonths=sum(numMonths),
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
  summarize(Tier=case_when(numYrs >= 5 & numWeeks >= 260 ~ "Level 2a",
                           numYrs >= 5 & numWeeks >= 130 ~ "Level 3a",
                           numYrs >= 5 & numMonths >= 48 ~ "Level 4a",
                           numYrs >= 5 & numQtrs >= 20 ~ "Level 5a",
                           # second set
                           numYrs >= 5 & numWeeks >= 260 ~ "Level 2b",
                           numYrs >= 5 & numWeeks >= 130 ~ "Level 3b",
                           numYrs >= 5 & numMonths >= 60 ~ "Level 4b",
                           numYrs >= 5 & numQtrs >= 20 ~ "Level 5b",
                           # third set
                           numYrs < 5 & numWeeks >= (numYrs*52) & max(lastYear)==2015 ~ "Level 2c",
                           numYrs < 5 & numWeeks >= (numYrs*26) & max(lastYear)==2015 ~ "Level 3c",
                           numYrs < 5 & numMonths >= (numYrs*12)  & max(lastYear)==2015 ~ "Level 4c",
                           numYrs < 5 & numQtrs >= (numYrs*4) & max(lastYear)==2015 ~ "Level 5c",
                           # fourth set
                           numYrs >= 3 & numWeeks >= 130 & numWeeks < 260 ~ "Level 2d",
                           numYrs >= 3 & numWeeks >= 65 & numWeeks < 130 ~ "Level 3d",
                           numYrs >= 3 & numMonths >= 24 & numMonths < 48 ~ "Level 4d",
                           numYrs >= 3 & numQtrs >= 10 & numQtrs < 20 ~ "Level 5d"))#,

Tiered_sites %>% group_by(Tier) %>% summarize(n=n())
# Write out results
write_feather(Tiered_sites,paste('4_gap_analysis/in/Tiers/',Con,'/',stateAb[1],'_Tiers.feather',sep=''))
write.csv(Tiered_sites,paste('4_gap_analysis/in/Tiers/',Con,'/',stateAb[1],'_Tiers.csv',sep=''),row.names=F,quote=F)

