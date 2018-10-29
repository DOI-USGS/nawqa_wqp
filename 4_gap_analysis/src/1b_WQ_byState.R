# GAP analysis script for the water quality coordination effort
# This script covers temperature, DO, pH, and SC

# To do - take out "DO saturation" con from the DO files

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

# Cons for this script - DO,Temp, pH, SC (conductivity)
Con<-"SC"
ConFull<-"conductivity"

# Crosswalk based on Agency name, NDA acroym, and WQX agency code - from Melissa
xWalk2<-loadWorkbook("4_gap_analysis/in/Organizations_WQP.xlsx") %>% readWorksheet(sheet="Organizations_WQP_testpull",header=T)

state2<-gsub('_',' ',state) # Does not work for DC
stateAb<-state.abb[grep(state2, state.name)][1]

# pull the files off of Google Drive - browser window for authentication should pop up the first time you do this
#indicator_files <- dir(paste0('1_wqpdata/out/data/',Con), pattern='*\\.feather.ind$', full.names=TRUE) %>%
#  grep(ConFull, ., value=TRUE)
#sapply(indicator_files, scipiper::gd_get, USE.NAMES=FALSE)

# read in the data and combine into one big file
#data_files <- sapply(indicator_files, scipiper::as_data_file, USE.NAMES=FALSE)
# Reading files locally
data_files <- dir(paste0('1_wqpdata/out/data/',Con), pattern=paste(state,'*', sep=''),full.names = T) %>%
 grep(paste('(',ConFull,')',sep=""), ., value=TRUE)
data_files<-data_files[!grepl('siteinfo|ind',data_files)] # can't do negative filters in lists with dplyr

# combine N/P data into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    # By State
    dat <- feather::read_feather(dfile) %>%
      # Below here are data coercion fixes for columns(1 - FL, 2- AL, 3- MA, 4 - MI, 5 - PA,6 - ?, 7-CO, 8/9-FL, 11-UT)
      mutate_at(vars(starts_with('ResultDepthHeightMeasure.MeasureValue'),starts_with('ActivityTopDepthHeightMeasure.MeasureValue'),
                     starts_with('ActivityBottomDepthHeightMeasure.MeasureValue'),starts_with('ResultAnalyticalMethod.MethodIdentifier')),
                funs(as.character)) %>%
      mutate(ConstituentGroup=ConFull)
  })
) %>% # Non-coercion changes to data, filling in NA, modifying some attributes
   # Drop date fields if NA
   drop_na(ActivityStartDate)

print (dim(combined_data))

conList<-combined_data %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName)%>%
  summarize(sum=length(OrganizationIdentifier))

# Looks like we need to add/address these orgs to xWalk2
print (conList$OrganizationIdentifier[!conList$OrganizationIdentifier %in% xWalk2$OrganizationIdentifier])
# Add the subsource ID to the table
conList<-left_join(conList,xWalk2,by='OrganizationIdentifier') %>%
  select(-one_of('OrganizationFormalName.y','natqw_source_state','natqw_source_type','source_name','subsource_name','natqw_sourceid','comment','Andy'))
write.table(conList,paste('4_gap_analysis/in/Org/',Con,'/',stateAb,'_OrgCon_List.csv',sep=''),row.names=F,sep=",") # changed to write.table at DC

# writeout data per State per Con (may not be necessary)
# The positive is the data type issues will be fixed, and data types consistent across fields
customFun  = function(DF) {
  #print(DF$ConstituentGroup)
  write.csv(DF,paste0("4_gap_analysis/in/Cons/",Con,"/",stateAb,"_",Con,".csv"))
  feather::write_feather(DF,paste("4_gap_analysis/in/Cons/",Con,"/",stateAb,"_",Con,".feather",sep=""))
  return(DF)
}
combined_data %>% group_by(ConstituentGroup) %>% do(customFun(.))

## Retrieve siteInfo data and built a siteFile
site_files <- dir(paste0('1_wqpdata/out/data/',Con), pattern=paste(state,'*',sep=""), full.names=TRUE) %>%
  grep('siteinfo', ., value=TRUE)
site_files<-site_files[!grepl('ind',site_files)] # can't do negative filters in lists with dplyr

# combine siteInfo data into one data frame and subset to unique rows
combined_sites <- bind_rows(
  lapply(site_files, function(dfile) {
    # By State
    cname <- gsub(paste0('1_wqpdata/out/data/',Con), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate_at(vars(starts_with('HorizontalAccuracyMeasure.MeasureValue'),starts_with('VerticalAccuracyMeasure.MeasureValue'),
                starts_with('ConstructionDateText'),starts_with('ContributingDrainageAreaMeasure.MeasureValue'),
                starts_with('VerticalMeasure.MeasureValue'),starts_with('VerticalMeasure.MeasureValue'),
                starts_with('SourceMapScaleNumeric')),funs(as.character))
  })
) %>% distinct()  %>% filter(dec_lon_va != 0 | dec_lat_va !=0) %>% # derive distinct rows, remove rows with no Lat/Lon data
  mutate(dec_lon_va=ifelse(dec_lon_va>0,-dec_lon_va,dec_lon_va)) # Re-calculate Lon as negative if it is positive

# write out site Info files
state<-gsub("_"," ",state)
stateAb<-state.abb[grep(state, state.name)][1]
feather::write_feather(combined_sites,paste("4_gap_analysis/in/Site/",Con,"/",stateAb,"_sites.feather",sep=""))
write.csv(combined_sites,paste("4_gap_analysis/in/Site/",Con,"/",stateAb,".csv",sep=""))

# write out information on missing sites (if applicable)
miss_staID<-combined_data$MonitoringLocationIdentifier[!combined_data$MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier]
write.csv(miss_staID,paste("4_gap_analysis/in/Site/",Con,"/",stateAb,"_MissingSites.csv",sep=""))

