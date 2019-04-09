# GAP analysis script for the water quality coordination effort
# This script covers Nutrients, and designtates fractions and
#      generic names

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
  state<-'Alabama'
}

# Cons for this script - Ammonia, nitrate, phosphorus, nitrogen, TKN
Con<-'TKN'
StartDate<-'1986-01-01'
EndDate<-'2015-12-31'

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

# If multiple species generalizations are present (i.e. Nitrate)
# Each list represents WQP characteristic names that fall under the generic Con/ConFull classification
if (Con=='Ammonia'){
  ConList<-c('Ammonia','Ammonia and ammonium','Ammonia as N','Ammonia as NH3','Ammonia as NH4','Ammonia uptake',
             'Ammonia-nitrogen','Ammonia-nitrogen as N','Ammonia-nitrogen as NH3','Ammonia-nitrogen as NH3',
             'Ammonium','Ammonium as N','Ammonium as NH4','Nitrogen, ammonia (NH3) + ammonium (NH4)',
             'Nitrogen, Ammonium (NH4) as N','Nitrogen, ammonium (NH4) as NH4')}

if (Con=='nitrate'){
  ConList<-c('Nitrate','Nitrite plus nitrate','Inorganic nitrogen (nitrate and nitrite)', # missing on first run
             'Inorganic nitrogen (nitrate and nitrite) as N','Nitrate + Nitrite','Nitrate as N','Nitrate as NO2',
             'Nitrate as NO3','Nitrate-N','Nitrite as N','Nitrogen, Nitrate (NO3) as NO3',
             'Nitrite','Nitrogen, Nitrite (NO2) as NO2','Nitrite as NO2')}

if (Con=='nitrogen'){
  ConList<-c('Inorganic nitrogen (ammonia, nitrate and nitrite','Inorganic nitrogen (NO2, NO3, & NH3)','Nitrogen',
             'Nitrogen as N','Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)','Nutrient-nitrogen',
             'Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3) as N','Nutrient-nitrogen as N',
             'Organic Nitrogen','Organic nitrogen','Organic Nitrogen as N','Particulate Organic Nitrogen and Particulate Nitrogen',
             'Total Nitrogen, mixed forms','Total Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)',
             'Total Particulate Nitrogen')}

if (Con=='phosphorus'){
  ConList<-c('Hydrolyzable phosphorus','Organic phosphorus','Organic phosphorus as P','Orthophosphate','Orthophosphate as P',
             'Orthophosphate as PO4','Polyphosphate','Soluble Reactive Phosphorus (SRP)','Soluble Reactive Phosphorus (SRP) as P',
             'Phosphate','Phosphate-phosphorus','Phosphate-phosphorus as P','Phosphate-phosphorus as PO4','Phosphorus',
             'Phosphorus as P','Phosphorus as PO4','Phosphorus, Particulate Organic','Phosphorus, total',
             'Total Phosphorus,  mixed forms','Particulate Inorganic Phosphorus',
             'Total Phosphorus, mixed forms',NA)}

if (Con=='TKN'){
  ConList<-c('Kjeldahl nitrogen','Kjeldahl nitrogen as N','Nitrogen, Ammonia + Organic','Nitrogen, Kjeldahl',
             'Total Kjeldahl nitrogen','Total Kjeldahl nitrogen (Organic N & NH3)')}

# Crosswalk based on Harmonization table
xWalk1<-loadWorkbook('4_gap_analysis/in/NH2012Trends_N.xlsx') %>% readWorksheet(sheet='data_harmonization',header=T) %>%
  filter(orig_parameter_name %in% ConList,qwdata_source %in% c('STORET/WQX','WQP'))

# Crosswalk based on Agency name, NDA acroym, and WQX agency code
xWalk2<-loadWorkbook('4_gap_analysis/in/Organizations_WQP.xlsx') %>% readWorksheet(sheet='WQP_Organizations_mlr',header=T)

#******************************
# Reading files locally
all_files <- dir('1_wqpdata/out/byState', pattern=paste0("^",state),full.names=T) %>% #pattern=paste0('*',state,'*'),full.names = T) %>%
  grep(paste('(',Con,')',sep=''), ., value=TRUE)
# Add ortho P data to the phosphorus files
if(Con=='phosphorus'){
  orthoPfiles<-dir('1_wqpdata/out/byState', pattern=paste0(state,'*'),full.names = T) %>%
    grep('(orthophosphate)', ., value=TRUE)
  all_files<-append(all_files,orthoPfiles)
}
data_files<-all_files[!grepl('siteinfo|ind',all_files)] # just the data files
site_files<-all_files[grepl('siteinfo|ind',all_files)] # just the site files

# combine into one data frame
combined_data <- bind_rows(
  lapply(data_files, function(dfile) {
    # By State
    dat <- feather::read_feather(dfile) %>%
      # Below here are data coercion fixes for columns
      mutate_at(vars(starts_with('ResultDepthHeightMeasure.MeasureValue'),starts_with('ActivityTopDepthHeightMeasure.MeasureValue'),
                     starts_with('ActivityBottomDepthHeightMeasure.MeasureValue'),starts_with('ResultAnalyticalMethod.MethodIdentifier'),
                     starts_with('ProjectIdentifier'),starts_with('SampleCollectionMethod.MethodIdentifier'),
                     starts_with('ActivityCommentText'),starts_with('LaboratoryName'),starts_with('ActivityDepthAltitudeReferencePointText')),
                list(as.character)) %>%
      mutate(ConstituentGroup=Con) %>%
      #See Read me for explanations on each constituent
      filter(!grepl('Quality Control',ActivityTypeCode) & ActivityMediaName %in% c('Water','water','Sediment') &
                ActivityMediaSubdivisionName %in% c('Surface Water',NA,'','Unknown') &
                !grepl('%',ResultMeasure.MeasureUnitCode) &
                ActivityStartDate > as.Date(StartDate) & ActivityStartDate < as.Date(EndDate))
  })
) %>%
  # Drop date fields if NA
  filter(!is.na(ActivityStartDate)) %>%
  mutate(ResultSampleFractionText=replace(ResultSampleFractionText,ResultSampleFractionText==" ",NA))


# subset by rows with data values in one of two fields
combined_data<-subset(combined_data,(!is.na(combined_data$ResultMeasureValue)) |
                        (!is.na(combined_data$DetectionQuantitationLimitMeasure.MeasureValue)))

print(combined_data %>% count(CharacteristicName))
print(combined_data %>% count(ResultMeasure.MeasureUnitCode))
print(combined_data %>% count(ResultSampleFractionText))
print (dim(combined_data))

#****************************
# Step 2 - nutrient crosswalk
# conList is a short table that is a unique list of the WQX agency code, Agency names, and Fraction
conList_MR<-combined_data  %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName,ResultSampleFractionText,
           ResultMeasure.MeasureUnitCode) %>%
  summarize(sum=length(OrganizationIdentifier))
# Org Lists and counts Requested by Melissa and Gretchen
write.csv(conList_MR,paste('4_gap_analysis/in/Org/',Con,'/',stateAb,'_OrgCon_List_MR.csv',sep=''),row.names=F,quote=F)

conList<-combined_data %>% mutate(OrganizationFormalName=gsub(',',' ',OrganizationFormalName)) %>%
  group_by(OrganizationIdentifier,OrganizationFormalName,CharacteristicName,ResultSampleFractionText)%>%
  summarize(sum=length(OrganizationIdentifier))
# print if any organizations not in Melissa's sheet
print (conList$OrganizationIdentifier[!conList$OrganizationIdentifier %in% xWalk2$OrganizationIdentifier])

orgs<-xWalk2[xWalk2$OrganizationIdentifier %in% conList$OrganizationIdentifier,]
# Add the subsource ID to the table
conList<-left_join(conList,orgs,by='OrganizationIdentifier') %>%
  select(-one_of('OrganizationFormalName.y','natqw_source_state','natqw_source_type','source_name','subsource_name','natqw_sourceid','comment','Andy'))
write.table(conList,paste('4_gap_analysis/in/Org/',Con,'/',stateAb,'_OrgCon_List.csv',sep=''),row.names=F,sep=',') # changed to write.table at DC

# Join constituents to Full agency information
# There was strange behavior for a few states when both join statements were merged
#       into a single pipeline
data_split_P<-combined_data %>% inner_join(xWalk2,by=c('OrganizationIdentifier')) %>% # Join to get the WQP Org Identifier
  inner_join(xWalk1,by=c('natqw_subsourceid'='natqw_subsourceid',
                          'CharacteristicName'='orig_parameter_name',
                          'ResultSampleFractionText'='orig_fraction'))
# LA DEQ data is counted twice for two sub-source ids linked to LA DEQ
if (stateAb=='LA'){
  data_split_P<-data_split_P %>% filter(natqw_subsourceid != 'LA DEQ WPD')
}

data_split_P %>% group_by(frac_mod) %>% summarize(count=n())

# original data
combined_data %>% group_by(OrganizationIdentifier) %>% summarize(count=n())
# Final join data
data_split_P %>% group_by(OrganizationIdentifier) %>% summarize(count=n())
# Count of different fraction methods
data_split_P %>% group_by(frac_mod) %>% summarize(count=n())

# combine siteInfo data into one data frame and subset to unique rows
combined_sites <- bind_rows(
  lapply(site_files, function(dfile) {
    # By State
    cname <- gsub(paste('1_wqpdata/out/data/',Con,'/',state,'_',sep=''), '', gsub('_00.\\.feather', '', dfile))
    dat <- feather::read_feather(dfile) %>%
      mutate_at(vars(starts_with('HorizontalAccuracyMeasure.MeasureValue'),starts_with('ConstructionDateText'),
                starts_with('VerticalMeasure.MeasureValue'),starts_with('VerticalAccuracyMeasure.MeasureValue'),
                starts_with('StateCode'),starts_with('SourceMapScaleNumeric'),
                starts_with('ContributingDrainageAreaMeasure.MeasureValue')),
                list(as.character)) %>%
      filter(!MonitoringLocationTypeName %in% c('Pond-Stormwater','Pond-Anchialine','Canal Transport','Canal Irrigation'))
  })
) %>% distinct()  %>% filter(LongitudeMeasure != 0 | LatitudeMeasure !=0) %>% # derive distinct rows, remove rows with no Lat/Lon data
  mutate(LongitudeMeasure=ifelse(LongitudeMeasure>0,-LongitudeMeasure,LongitudeMeasure)) # Re-calculate Lon as negative if it is positive

# Ensure sites are same between both datasets
data_split_P<-data_split_P %>% filter(MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier)

# The positive is the data type issues will be fixed, and data types consistent across fields
customFun  = function(DF) {
  write.csv(DF,paste0('4_gap_analysis/in/Cons/',Con,'/',stateAb,'_',unique(DF$ConstituentGroup),'.csv'))
  feather::write_feather(DF,paste('4_gap_analysis/in/Cons/',Con,'/',stateAb,'_',unique(DF$ConstituentGroup),'.feather',sep=''))
  return(DF)
}
data_split_P %>% group_by(ConstituentGroup) %>% do(customFun(.))

combined_sites<-combined_sites %>% filter(MonitoringLocationIdentifier %in% combined_data$MonitoringLocationIdentifier)

feather::write_feather(combined_sites,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'_sites.feather',sep=''))
write.csv(combined_sites,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'.csv',sep=''))

# write out information on missing sites (if applicable)
miss_staID<-combined_data$MonitoringLocationIdentifier[!combined_data$MonitoringLocationIdentifier %in% combined_sites$MonitoringLocationIdentifier]
if (length(miss_staID) > 0){
  write.csv(miss_staID,paste('4_gap_analysis/in/Site/',Con,'/',stateAb,'_MissingSites.csv',sep=''))
}


