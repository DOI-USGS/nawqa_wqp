# This file contains a first cut at partitioning whatever data are currently in
# your 1_wqpdata/out/data folder into these scenarios:
#   A (ideal) - sites with 12 samples per year, dissolved P, TP, and NO3
#   B (almost there, need more constituents) - sites with 12 samples per year, 1 or two of dissolved P, TP, NO3
#   C (almost there, need more frequency) - sites with 6 samples per year, dissolved P, TP, NO3
# (as you'll see, I've had to back the number of constituents down to 2 for A
# and C to get any sites for AL)

library(tidyverse)
library(feather) # install with install.packages('feather')
library(scipiper) # install with devtools::install_github('USGS-R/scipiper')

state<-"Florida"

# pull the files off of Google Drive - browser window for authentication should
# pop up the first time you do this
indicator_files <- dir('1_wqpdata/out/data', pattern='*\\.feather.ind$', full.names=TRUE) %>%
  grep('(nitrate|phosphorus)', ., value=TRUE)
sapply(indicator_files, scipiper::gd_get, USE.NAMES=FALSE)

# read in the data and combine into one big file
#data_files <- sapply(indicator_files, scipiper::as_data_file, USE.NAMES=FALSE)
data_files <- dir('1_wqpdata/out/data', pattern=paste(state,'*', sep=''),full.names = T) %>%
  grep('(nitrate|phosphorus)', ., value=TRUE)

# combined data may not work
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    cname <- gsub(paste('1_wqpdata/out/data/',state,'_',sep=""), '', gsub("_00.\\.feather", '', dfile))
    dat <- feather::read_feather(dfile) %>%
     # mutate(ConstituentGroup=cname) %>%
    # Below here are data conversion problems
      mutate_at(vars(starts_with("Result")),funs(as.character)) %>% # FL
      mutate_at(vars(ends_with("HeightMeasure.MeasureValue")),funs(as.character))# AL
  })
)%>%
  mutate(ResultSampleFractionText=ifelse(is.na(ResultSampleFractionText),"Unknown",ResultSampleFractionText))

# If we are using the ActivityStartDateTime field, get rid of NAs
# Don't understand the replacement error warning error, lens of both
# objects are 1
combined_data$ActivityStartDateTime[is.na(combined_data$ActivityStartDateTime)]<-
  as.POSIXlt(paste(combined_data$ActivityStartDate[is.na(combined_data$ActivityStartDateTime)],
                   "12:00:00",sep=" "),format="%Y-%m-%d %H:%M:%S")

# Crosswalk implementation, check with Gretchen
# Does not include Organic P (236) and Polyphosphate (55)
# List of Agency Names, Constitutents, and Filter Type from WQ Pull
combined_data$OrganizationIdentifier[grepl("USGS*", combined_data$OrganizationIdentifier)]<-"USGS"
# Subset conList to only P
conList<-filter(combined_data,grepl("*(P|p)hos*",CharacteristicName)) %>% # Check for Orthophosphate and SRP with this
  distinct(OrganizationIdentifier,CharacteristicName,ResultSampleFractionText) %>% #6/28 removed organization Name OrganizationFormalName,
mutate(ResultSampleFractionText=ifelse(is.na(ResultSampleFractionText),"Unknown",ResultSampleFractionText))

# Crosswalk based on Harmonization table
xWalk1<-read.csv("andy/NH2012.csv",stringsAsFactors = F)
# Crosswalk based on Agency List and acronyms
xWalk2<-read.csv("andy/Crosswalk.csv",stringsAsFactors = F)

# Join constituents to Full agency information
conAgency<-inner_join(x=conList,y=xWalk2,by=c("OrganizationIdentifier"="WQX"))
conSampInfo<-inner_join(x=conAgency,y=xWalk1,by=c("SIR"="natqw_subsourceid",
                                                  "CharacteristicName"="orig_parameter_name",
                                                  "ResultSampleFractionText"="orig_fraction")) %>%
  distinct(OrganizationIdentifier,CharacteristicName,ResultSampleFractionText,Final)#added Fraction 6/27

data_split_P<-combined_data %>% left_join(x=combined_data,y=conSampInfo,
                                          by=c("OrganizationIdentifier","CharacteristicName",
                                               "ResultSampleFractionText")) %>%
  mutate(Final=replace(Final,is.na(Final),"Nitrate")) %>% rename(ConstituentGroup=Final)

data_split_P %>% group_by(ConstituentGroup) %>% summarise(NumValues=length(ConstituentGroup))

#blah<-blah %>% select(OrganizationIdentifier,OrganizationFormalName,Final,CharacteristicName,ResultSampleFractionText)

UK<-data_split_P[data_split_P$ResultSampleFractionText=="Unknown",]
UK %>% group_by(ConstituentGroup) %>% summarise (NumValues=length(ConstituentGroup))

#hoho<-data_split_P %>% group_by(MonitoringLocationIdentifier,CharacteristicName) %>% summarize(num=length(CharacteristicName))

# count the number of target constituents (should be TP, dissolved P, and NO3,
# but for now we'll just do any-P and NO3) per sample. this takes a few seconds,
# watch the date format, it will not mutate if the date is not a POSIX
constituents_per_sample <- data_split_P %>%
  select(ConstituentGroup,MonitoringLocationIdentifier,ResultMeasureValue,ActivityStartDateTime) %>%
  mutate(Year = as.integer(format(ActivityStartDateTime, '%Y'))) %>%
  group_by(MonitoringLocationIdentifier, ActivityStartDateTime, Year) %>%
  summarize(NumConstituentsPerSample=length(unique(ConstituentGroup))) %>%
  ungroup()

# count the number of samples per year, separating into samples with one
# constituent or two
samples_per_year <- constituents_per_sample %>%
  group_by(MonitoringLocationIdentifier, Year) %>%
  summarize(
    NumSamplesPerYear_Min1Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 1]),
    NumSamplesPerYear_Min2Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 2]),
    NumSamplesPerYear_Min3Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 3])) %>%
  ungroup() %>%
  gather(temporary, NumSamplesPerYear, starts_with('NumSamplesPerYear')) %>%
  extract(temporary, 'MinNumConstituents', "NumSamplesPerYear_Min([[:alnum:]])Const") %>%
  mutate(MinNumConstituents=as.integer(MinNumConstituents))

# compute sampling frequencies for samples containing at least 1, at least 2, or
# 3 different constitutent types
sample_rates <- samples_per_year %>%
  # make sure all 6 years are represented for each site, even if there were 0
  # samples some years
  filter(Year %in% 2010:2015) %>%
  complete(Year = 2010:2015) %>%
  mutate(NumSamplesPerYear = ifelse(is.na(NumSamplesPerYear), 0, NumSamplesPerYear)) %>%
  # compute various stats for frequency over all target years
  group_by(MonitoringLocationIdentifier, MinNumConstituents) %>%
  summarize(
    MeanSamplesPerYear=mean(NumSamplesPerYear),
    MinSamplesPerYear=min(NumSamplesPerYear),
    MaxSamplesPerYear=max(NumSamplesPerYear),
    MedianSamplesPerYear=median(NumSamplesPerYear))

# split into scenarios
# Do some more work here
scenario_empty <- sample_rates %>%
  filter(MinNumConstituents == 3, MinSamplesPerYear >= 12)
scenario_A <- sample_rates %>%
  filter(MinNumConstituents == 2, MinSamplesPerYear >= 12)
scenario_B <- sample_rates %>%
  filter(MinNumConstituents == 1, MinSamplesPerYear >= 12)
scenario_C <- sample_rates %>%
  filter(MinNumConstituents == 2, MinSamplesPerYear >= 6)

siteInfo<-rbind(scenario_A,scenario_B,scenario_C)
siteInfo$cat<-ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_A$MonitoringLocationIdentifier,"A) 2 cons., >=12 samples/year",
                     ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_B$MonitoringLocationIdentifier,"B) 1 cons., >=12 samples/year",
                            ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_C$MonitoringLocationIdentifier,"C) 2 cons. >=6 samples/year","D) <6 samples/year")))

state<-gsub("_"," ",state)
stateAb<-state.abb[grep(state, state.name)]

write_feather(siteInfo,paste("andy/",stateAb,"_samples.feather",sep=""))
write.csv(siteInfo,paste("andy/",stateAb,"_samples.csv",sep=""),row.names=F,quote=F)

# #### explore results ####
#
# # there are at least a few thousand of each, many fewer for dissolved P
# #data_split_P$ConstituentGroup %>% table
#
# # site counts for each scenario - A and B are relatively infrequent, but C is
# # pretty common
# nrow(scenario_A)
# nrow(scenario_B)
# nrow(scenario_C)
#
# # histogram of samples per year grouped by minimum number of constituents per
# # sample - shows that samples with 3 constituents don't currently exist in AL,
# # samples that have 1 almost always have 2 constituents, and many more sites
# # have 5-10 samples per year than >= 12
# sample_rates %>%
#   mutate(MinNumConstituents=as.factor(MinNumConstituents)) %>%
#   ggplot(aes(x=MeanSamplesPerYear, group=MinNumConstituents, fill=MinNumConstituents)) +
#   geom_histogram(binwidth=1, position='dodge')


