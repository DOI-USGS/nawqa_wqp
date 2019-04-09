# GAP analysis script for the water quality coordination effort
# This script performs the initial processing and filtering of data

library(tidyverse)
library(feather) # install with install.packages('feather')
library(scipiper) # install with devtools::install_github('USGS-R/scipiper')
library(XLConnect) # read xlsx files

#***********************
# batch mode or manually
runMode<-'batch'
if (runMode=='batch'){
  setwd('D:/abock/NAWQA/WQR/nawqa_wqp') # Set locally for running in batch mode for quick processing of all states
  args<-commandArgs(trailingOnly = T)
  state<-args[1]
}else{
  state<-'Kentucky'
}

# Constituent this script runs for
Con<-'DO'
StartDate<-'1986-01-01'
EndDate<-'2015-12-31'

# Crosswalk based on Agency name, NDA acroym, and WQX agency code
xWalk2<-loadWorkbook('4_gap_analysis/in/Organizations_WQP.xlsx') %>% readWorksheet(sheet='WQP_Organizations_mlr',header=T)

# Designate state abbreviation
state2<-gsub('_',' ',state)
stateAb<-state.abb[grep(state2, state.name)][1]
# State abbreviations for states not covered by R state abbreviation database
if(state == 'Northern_Mariana_Islands') stateAb='MP'
if(state == 'Guam') stateAb='GU'
if(state == 'American_Samoa') stateAb='AS'
if(state == 'District_of_Columbia') stateAb='DC'
if(state == 'Puerto_Rico') stateAb='PR'
if(state == 'Virgin_Islands') stateAb='VI'

# Reading files locally
all_files <- dir('1_wqpdata/out/byState', pattern=paste0("^",state),full.names=T) %>% #pattern=paste0('*',state,'*'),full.names = T) %>%
 grep(paste('(',Con,')',sep=''), ., value=TRUE)
data_files<-all_files[!grepl('siteinfo|ind',all_files)] # just the data files
site_files<-all_files[grepl('siteinfo|ind',all_files)] # just the site files

# combine files into one data frame for each state/CON combination
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    dat <- feather::read_feather(dfile) %>%
      # Below here are data coercion fixes for the specific columns
      mutate_at(vars(starts_with('ResultDepthHeightMeasure.MeasureValue'),starts_with('ActivityTopDepthHeightMeasure.MeasureValue'),
                     starts_with('ActivityBottomDepthHeightMeasure.MeasureValue'),starts_with('ResultAnalyticalMethod.MethodIdentifier'),
                     starts_with('ProjectIdentifier'),starts_with('SampleCollectionMethod.MethodIdentifier'),
                     starts_with('ActivityCommentText'),starts_with('LaboratoryName')),
                list(as.character)) %>%
      mutate(ConstituentGroup=Con) %>%
      # filter out quality control rows and media types
      filter(!grepl('Quality Control',ActivityTypeCode) & ActivityMediaName %in% c('Water','water','Sediment') &
              ActivityMediaSubdivisionName %in% c('Surface Water',NA,'','Unknown') &
              # Filter out DO % saturation, not confident we can use this information without pressure and temperature
              !ResultMeasure.MeasureUnitCode %in% c('%','% saturatn','%         ','deg C',' ','None','count','m','nu') &
               !is.na(ResultMeasure.MeasureUnitCode)&
               !CharacteristicName %in% c('Dissolved oxygen saturation') &
              ActivityStartDate > as.Date(StartDate) & ActivityStartDate < as.Date(EndDate))
  })
) %>% # Activity Start Date is not completed, can't be confident when measurement was actually taken
   filter(!is.na(ActivityStartDate))

# Both these fields may have valid measurement or concentration values
# This keeps rows that have data in or the other
combined_data<-subset(combined_data,(!is.na(combined_data$ResultMeasureValue)) |
                        (!is.na(combined_data$DetectionQuantitationLimitMeasure.MeasureValue)))

# Name and unit summaries
print(combined_data %>% count(CharacteristicName))
print(combined_data %>% count(ResultMeasure.MeasureUnitCode))
print (dim(combined_data))

# Create organization list and
conList<-combined_data %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName)%>%
  summarize(sum=length(OrganizationIdentifier))
# Orgs to add/address these orgs to xWalk2
print (conList$OrganizationIdentifier[!conList$OrganizationIdentifier %in% xWalk2$OrganizationIdentifier])

# Add the subsource ID to the table, remove duplicate columns
conList<-left_join(conList,xWalk2,by='OrganizationIdentifier') %>%
  select(-one_of('OrganizationFormalName.y','natqw_source_state','natqw_source_type','source_name','subsource_name','natqw_sourceid','comment','Andy'))
write.table(conList,paste('4_gap_analysis/in/Org/',Con,'/',stateAb,'_OrgCon_List.csv',sep=''),row.names=F,sep=',') # changed to write.table at DC

# combine siteInfo data into one data frame and subset to unique rows
combined_sites <- bind_rows(
  lapply(site_files, function(dfile) {
    # By State
    cname <- gsub(paste0('1_wqpdata/out/data/',Con), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate_at(vars(starts_with('HorizontalAccuracyMeasure.MeasureValue'),starts_with('VerticalAccuracyMeasure.MeasureValue'),
                starts_with('ConstructionDateText'),starts_with('ContributingDrainageAreaMeasure.MeasureValue'),
                starts_with('VerticalMeasure.MeasureValue'),starts_with('VerticalMeasure.MeasureValue'),
                starts_with('SourceMapScaleNumeric'),starts_with('DrainageAreaMeasure.MeasureValue'),
                starts_with('StateCode')),list(as.character)) %>%
      filter(!MonitoringLocationTypeName %in% c('Pond-Stormwater','Pond-Anchialine','Canal Transport','Canal Irrigation'))
  })
) %>% distinct()  %>% filter(LongitudeMeasure != 0 | LatitudeMeasure !=0) %>% # derive distinct rows, remove rows with no Lat/Lon data
  mutate(LongitudeMeasure=ifelse(LongitudeMeasure>0,-LongitudeMeasure,LongitudeMeasure)) # Re-calculate Lon as negative if it is positive


# Ensure sites are same between both datasets
combined_data<-combined_data %>% filter(MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier)
customFun  = function(DF) {
  write.csv(DF,paste0('4_gap_analysis/in/Cons/',Con,'/',stateAb,'_',Con,'.csv'))
  feather::write_feather(DF,paste('4_gap_analysis/in/Cons/',Con,'/',stateAb,'_',Con,'.feather',sep=''))
  return(DF)
}
combined_data %>% group_by(ConstituentGroup) %>% do(customFun(.))

combined_sites<-combined_sites %>% filter(MonitoringLocationIdentifier %in% combined_data$MonitoringLocationIdentifier)
# write out site Info files
feather::write_feather(combined_sites,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'_sites.feather',sep=''))
write.csv(combined_sites,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'.csv',sep=''))

# write out information on missing sites (if applicable)
miss_staID<-combined_data$MonitoringLocationIdentifier[!combined_data$MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier]
if (length(miss_staID) > 0){
  write.csv(miss_staID,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'_MissingSites.csv',sep=''))
}
