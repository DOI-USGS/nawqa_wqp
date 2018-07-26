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

state<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state, state.name)][1]

# Crosswalk based on Harmonization table
xWalk1<-loadWorkbook("4_gap_analysis/in/NH2012Trends.xlsx") %>% readWorksheet(sheet="data_harmonization",header=T)
# Crosswalk based on Agency name, NDA acroym, and WQX agency code - from Melissa
xWalk2<-loadWorkbook("4_gap_analysis/in/Organizations_WQP_testPull.xlsx") %>% readWorksheet(sheet="Organizations_WQP_testpull",header=T)

# pull the files off of Google Drive - browser window for authentication should pop up the first time you do this
indicator_files <- dir('1_wqpdata/out/data', pattern='*\\.feather.ind$', full.names=TRUE) %>%
  grep('(nitrate|phosphorus)', ., value=TRUE)
sapply(indicator_files, scipiper::gd_get, USE.NAMES=FALSE)

# read in the data and combine into one big file
data_files <- sapply(indicator_files, scipiper::as_data_file, USE.NAMES=FALSE)
# Reading files locally
#data_files <- dir('1_wqpdata/out/data', pattern=paste(state,'*', sep=''),full.names = T) %>%
#  grep('(nitrate|phosphorus)', ., value=TRUE)

# combine N/P data into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    # By State
    cname <- gsub(paste('1_wqpdata/out/data/',state,'_',sep=''), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate(ConstituentGroup=cname) %>%
      # Below here are data coercion fixes for columns(1 - FL, 2- AL, 3- MA, 4 - MI, 5 - PA )
      mutate_at(vars(starts_with('Result'),ends_with('HeightMeasure.MeasureValue'),ends_with('HeightMeasure.MeasureValue'),
                     starts_with('Measure'),starts_with('LaboratoryName'),starts_with('ActivityComment')),
                funs(as.character))
  })
) %>% # Non-coercion changes to data, filling in NA, modifying some attributes
   mutate(ResultSampleFractionText=ifelse(is.na(ResultSampleFractionText),'Unknown',ResultSampleFractionText)) %>%
   # Remove nitrite samples, likely groundwater
   filter(!CharacteristicName == 'Nitrite') %>%
   # Drop date fields if NA
   drop_na(ActivityStartDate)

#****************************
# Step 2 - nutrient crosswalk
# conList is a short table that is a unique list of the WQX agency code, Agency names, and Fraction
conList<-combined_data  %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName,ResultSampleFractionText) %>%
  summarize(sum=length(OrganizationIdentifier))
# Requested by Melissa and Gretchen
write.csv(conList,paste('andy/Org/',stateAb,'_OrgCon_List.csv',sep=''),row.names=F,quote=F)

# Join constituents to Full agency information
data_split_P<-combined_data %>% filter(grepl('*(P|p)hos*',CharacteristicName)) %>%
  inner_join(xWalk2,by=c('OrganizationIdentifier')) %>% # Join to get the WQP Org Identifier
  inner_join(xWalk1,by=c('natqw_subsourceid'='natqw_subsourceid',
                          'CharacteristicName'='orig_parameter_name',
                          'ResultSampleFractionText'='orig_fraction'))

#Printed data summaries, number of values per constituent group
data_split_P %>% group_by(ConstituentGroup) %>% summarise(NumValues=length(ConstituentGroup))
# Number of instances where the filtering method is 'unknown'
UK<-data_split_P[data_split_P$ResultSampleFractionText=='Unknown',]
UK %>% group_by(ConstituentGroup) %>% summarise (NumValues=length(ConstituentGroup))

# writeout data per State per Con (may not be necessary)
# The positive is the data type issues will be fixed, and data types consistent across fields
customFun  = function(DF) {
  print(DF$ConstituentGroup)
  write.csv(DF,paste0("andy/Cons/",stateAb,"_",unique(DF$ConstituentGroup),".csv"))
  return(DF)
}
combined_data %>% group_by(ConstituentGroup) %>% do(customFun(.))

#****************************
# Step 3 - Filtering methods

# count the number of target constituents per sample. this takes a few seconds,
# watch the date format, it will not mutate if the date is not a POSIX
constituents_per_sample <- data_split_P %>%
  select(ConstituentGroup,MonitoringLocationIdentifier,ResultMeasureValue,ActivityStartDate) %>%
  mutate(Year = as.integer(format(ActivityStartDate, '%Y'))) %>%
  group_by(MonitoringLocationIdentifier, ActivityStartDate, Year) %>%
  summarize(NumConstituentsPerSample=length(unique(ConstituentGroup))) %>%
  ungroup()

# count the number of samples per year, separating into samples with one constituent or two
samples_per_year <- constituents_per_sample %>%
  group_by(MonitoringLocationIdentifier, Year) %>%
  summarize(
    NumSamplesPerYear_Min1Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 1]),
    NumSamplesPerYear_Min2Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 2]),
    NumSamplesPerYear_Min3Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 3])) %>%
  ungroup() %>%
  gather(temporary, NumSamplesPerYear, starts_with('NumSamplesPerYear')) %>%
  extract(temporary, 'MinNumConstituents', 'NumSamplesPerYear_Min([[:alnum:]])Const') %>%
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
    MedianSamplesPerYear=median(NumSamplesPerYear),
    NoYears=n_distinct(Year))

# split into scenarios
# Do some more work here
scenario_A <- sample_rates %>%
  filter(MinNumConstituents == 3, MinSamplesPerYear >= 12, NoYears>4)
scenario_B <- sample_rates %>%
  filter(MinNumConstituents < 3, MinSamplesPerYear >= 12, NoYears>4)
scenario_C <- sample_rates %>%
  filter(MinNumConstituents == 3, MinSamplesPerYear >= 6 & MinSamplesPerYear < 12, NoYears>4)
scenario_D <- sample_rates %>%
  filter(MinNumConstituents < 3, MinSamplesPerYear >= 6 & MinSamplesPerYear < 12, NoYears>4)
scenario_E <- sample_rates %>%
  filter(MinNumConstituents == 3, MinSamplesPerYear >= 4 & MinSamplesPerYear < 6, NoYears>4)

siteInfo<-rbind(scenario_A,scenario_B,scenario_C,scenario_D,scenario_E)
siteInfo$cat<-ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_A$MonitoringLocationIdentifier,'A) 3 constituents (con.), >= 12 samples/year',
                     ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_B$MonitoringLocationIdentifier,'B) 1-2 con., >= 12 samples/year',
                            ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_C$MonitoringLocationIdentifier,'C) 3 con. 6-12 samples/year',
                                   ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_D$MonitoringLocationIdentifier,'D) 1-2 con. 6-12 samples/year',
                                          ifelse(siteInfo$MonitoringLocationIdentifier %in% scenario_E$MonitoringLocationIdentifier,'E) 3 con. 4-6 samples/year',NA)
                            ))))



write_feather(siteInfo,paste('andy/',stateAb,'_samples.feather',sep=''))
write.csv(siteInfo,paste('andy/',stateAb,'_samples.csv',sep=''),row.names=F,quote=F)

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


