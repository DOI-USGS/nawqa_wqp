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

state2<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state2, state.name)][1]

# Crosswalk based on Harmonization table
xWalk1<-loadWorkbook("4_gap_analysis/in/NH2012Trends.xlsx") %>% readWorksheet(sheet="data_harmonization",header=T)
# Crosswalk based on Agency name, NDA acroym, and WQX agency code - from Melissa
xWalk2<-loadWorkbook("4_gap_analysis/in/Organizations_WQP_testPull.xlsx") %>% readWorksheet(sheet="Organizations_WQP_testpull",header=T)

# # pull the files off of Google Drive - browser window for authentication should pop up the first time you do this
# indicator_files <- dir('1_wqpdata/out/data', pattern='*\\.feather.ind$', full.names=TRUE) %>%
#   grep('(nitrate|phosphorus)', ., value=TRUE)
# sapply(indicator_files, scipiper::gd_get, USE.NAMES=FALSE)

# read in the data and combine into one big file
#data_files <- sapply(indicator_files, scipiper::as_data_file, USE.NAMES=FALSE)
# Reading files locally
data_files <- dir('1_wqpdata/out/data', pattern=paste(state,'*', sep=''),full.names = T) %>%
 grep('(nitrate|phosphorus)', ., value=TRUE)
data_files<-data_files[!grepl('siteinfo',data_files)] # can't do negative filters in lists with dplyr

# combine N/P data into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    # By State
    cname <- gsub(paste('1_wqpdata/out/data/',state,'_',sep=''), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate(ConstituentGroup=cname) %>%
      # Below here are data coercion fixes for columns(1 - FL, 2- AL, 3- MA, 4 - MI, 5 - PA,6 - ?, 7-CO, 8/9-FL)
      mutate_at(vars(starts_with('Result'),ends_with('HeightMeasure.MeasureValue'),ends_with('HeightMeasure.MeasureValue'),
                     starts_with('Measure'),starts_with('LaboratoryName'),starts_with('ActivityComment'),
                     starts_with('HorizontalAccuracyMeasure.MeasureValue'),starts_with('ConstructionDateText'),
                     starts_with('VerticalMeasure.MeasureValue'),starts_with('VerticalAccuracyMeasure.MeasureValue')),
                funs(as.character))
  })
) %>% # Non-coercion changes to data, filling in NA, modifying some attributes
   mutate(ResultSampleFractionText=ifelse(is.na(ResultSampleFractionText),'Unknown',ResultSampleFractionText)) %>%
   # Remove nitrite samples, likely groundwater
   filter(!CharacteristicName == 'Nitrite') %>%
   # Drop date fields if NA
   drop_na(ActivityStartDate)

print (dim(combined_data))

#****************************
# Step 2 - nutrient crosswalk
# conList is a short table that is a unique list of the WQX agency code, Agency names, and Fraction
conList_MR<-combined_data  %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName,ResultSampleFractionText,
           ResultMeasure.MeasureUnitCode) %>%
  summarize(sum=length(OrganizationIdentifier))
# Requested by Melissa and Gretchen
write.csv(conList_MR,paste('4_gap_analysis/in/Org/',stateAb,'_OrgCon_List_MR.csv',sep=''),row.names=F,quote=F)

conList<-combined_data %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName,ResultSampleFractionText)%>%
  summarize(sum=length(OrganizationIdentifier))
# print if any organizations not in Melissa's sheet
print (conList$OrganizationIdentifier[!conList$OrganizationIdentifier %in% xWalk2$OrganizationIdentifier])
orgs<-xWalk2[xWalk2$OrganizationIdentifier %in% conList$OrganizationIdentifier,]
# Add the subsource ID to the table
conList<-left_join(conList,orgs,by='OrganizationIdentifier') %>%
  select(-one_of('OrganizationFormalName.y','natqw_source_state','natqw_source_type','source_name','subsource_name','natqw_sourceid','comment','Andy'))
write.table(conList,paste('4_gap_analysis/in/Org/',stateAb,'_OrgCon_List.csv',sep=''),row.names=F,sep=",") # changed to write.table at DC

# Join constituents to Full agency information
data_split_P<-combined_data %>% filter(grepl('*(P|p)hos*',CharacteristicName)) %>%
  inner_join(xWalk2,by=c('OrganizationIdentifier')) %>% # Join to get the WQP Org Identifier
  inner_join(xWalk1,by=c('natqw_subsourceid'='natqw_subsourceid',
                          'CharacteristicName'='orig_parameter_name',
                          'ResultSampleFractionText'='frac_mod'))

#Printed data summaries, number of values per constituent group
data_split_P %>% group_by(ConstituentGroup) %>% summarise(NumValues=length(ConstituentGroup))
# Number of instances where the filtering method is 'unknown'
UK<-data_split_P[data_split_P$ResultSampleFractionText=='Unknown',]
UK %>% group_by(ConstituentGroup) %>% summarise (NumValues=length(ConstituentGroup))

combined_data2<-bind_rows(combined_data[combined_data$ConstituentGroup=="nitrate",],data_split_P)

# writeout data per State per Con (may not be necessary)
# The positive is the data type issues will be fixed, and data types consistent across fields
customFun  = function(DF) {
  #print(DF$ConstituentGroup)
  write.csv(DF,paste0("4_gap_analysis/in/Cons/",stateAb,"_",unique(DF$ConstituentGroup),".csv"))
  feather::write_feather(DF,paste("4_gap_analysis/in/Cons/",stateAb,"_",unique(DF$ConstituentGroup),".feather",sep=""))
  return(DF)
}
combined_data2 %>% group_by(ConstituentGroup) %>% do(customFun(.))

## Retrieve siteInfo data and built a siteFile
site_files <- dir('1_wqpdata/out/data', pattern=paste(state,'*',sep=""), full.names=TRUE) %>%
  grep('siteinfo', ., value=TRUE)

# combine siteInfo data into one data frame and subset to unique rows
combined_sites <- bind_rows(
  lapply(site_files, function(dfile) {
    # By State
    cname <- gsub(paste('1_wqpdata/out/data/',state,'_',sep=''), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate_at(vars(starts_with('HorizontalAccuracyMeasure.MeasureValue'),starts_with('ConstructionDateText'),
                starts_with('VerticalMeasure.MeasureValue'),starts_with('VerticalAccuracyMeasure.MeasureValue'),
                starts_with('StateCode')),
                funs(as.character))
  })
) %>% distinct()  %>% filter(dec_lon_va != 0 | dec_lat_va !=0) %>% # derive distinct rows, remove rows with no Lat/Lon data
  mutate(dec_lon_va=ifelse(dec_lon_va>0,-dec_lon_va,dec_lon_va)) # Re-calculate Lon as negative if it is positive

print(class(combined_sites$StateCode))

# write out site Info files
state<-gsub("_"," ",state)
stateAb<-state.abb[grep(state, state.name)][1]
feather::write_feather(combined_sites,paste("4_gap_analysis/in/Site/",stateAb,"_sites.feather",sep=""))
write.csv(combined_sites,paste("4_gap_analysis/in/Site/",stateAb,".csv",sep=""))

# write out information on missing sites (if applicable)
miss_staID<-combined_data$MonitoringLocationIdentifier[!combined_data$MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier]
write.csv(miss_staID,paste("4_gap_analysis/in/Site/",stateAb,"_MissingSites.csv",sep=""))
