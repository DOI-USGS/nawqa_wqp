# This identifies closest streamgage and NHDflowline to each WQM site

library(tidyverse)
library(RANN)
library(lubridate)
library(zoo)

# Select the constituent Group
Con<-"Nut"

# Read in files
# NHD Directory, lists the regional directorys
nhdDir<-list.dirs("4_gap_analysis/in/GIS",recursive=F)
# Site Files by state
siteFiles<-dir(paste("4_gap_analysis/in/Site/",Con,sep=""), pattern="*sites.feather",full.names = T)
# Gages derived from an XML file
gages<-feather::read_feather(paste0("4_gap_analysis/in/Site/",Con,"/ClosestGage.feather"))

# for each gage derive flow stats
if (!file.exists(paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))){
  options(warn=1)
  FDC <- bind_rows(
    lapply(1:length(gages$STA), function(i) {
      gage<-gages$site_no[i]
      #print(i)
      if (file.exists(paste("4_gap_analysis/in/Q/bySTA/",gage,".csv",sep=""))){
        # read in Q file
        inQ<-read.csv(paste("4_gap_analysis/in/Q/bySTA/",gage,".csv",sep=""),colClasses="character")
        # Test for non-regular headings, convert if necessary
        if("X_.Final._00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_.Final._00060_00003}
        if("X_.Regression._00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_.Regression._00060_00003}
        if("X_.Primary.Stream.Flow._00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_.Primary.Stream.Flow._00060_00003}
        if("X_Red.River.Only_00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_Red.River.Only_00060_00003}
        if("X_FROM.DCP_00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_FROM.DCP_00060_00003}
        if("X_Powerhouse.Releases_00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_Powerhouse.Releases_00060_00003}
        if("X_.Data.from.10.1.1992.Forward._00060_00003" %in% colnames(inQ)){
          inQ$X_00060_00003<-inQ$X_.Data.from.10.1.1992.Forward._00060_00003}

        # Keep just minimal columns
        inQ<-inQ[c("agency_cd","site_no","dateTime","X_00060_00003")]
        inQ$dateTime<-as.Date(inQ$dateTime)
        inQ$X_00060_00003<-as.numeric(inQ$X_00060_00003)
        inQ$WaterYear<-ifelse(month(inQ$dateTime)<10,year(inQ$dateTime),year(inQ$dateTime)+1)
        inQ$MonthDay<-paste0(months(inQ$dateTime)," ",mday(inQ$dateTime))

        # get Q from each station for each date
        # 3 gages have NAs for entire time series and will generate soft warnings
        # stats
        # 1 - FDC (entire POR); 2-FDC (by year); 3 - 7-day MA; 4 - 15-day MA; 5 - 31-day MA;
        # 6 - 7-day rank; 7 - 15-day rank; 8 - 31-day rank
        nearQ_singleGage<-inQ %>% filter(site_no==as.character(gage)) %>% mutate(fdr=hydroTSM::fdc(X_00060_00003,plot=F)) %>%
          group_by(WaterYear) %>% mutate(fdrAnn=hydroTSM::fdc(X_00060_00003,plot=F)) %>% ungroup() %>%
          mutate(ma7=rollmean(X_00060_00003,7,malign='center',fill=NA),ma15=rollmean(X_00060_00003,15,malign='center',fill=NA),
                 ma31=rollmean(X_00060_00003,31,malign='center',fill=NA)) %>%
          mutate(rank7=zoo::rollapply(inQ$X_00060_00003,7,rank,fill=NA,align="center")[,4],
                 rank15=zoo::rollapply(inQ$X_00060_00003,15,rank,fill=NA,align="center")[,8],
                 rank31=zoo::rollapply(inQ$X_00060_00003,31,rank,fill=NA,align="center")[,16]) %>%
          group_by(MonthDay) %>%
          mutate(dayMedian=X_00060_00003/median(X_00060_00003,na.rm=T),
                 dayMean=X_00060_00003/mean(X_00060_00003,na.rm=T))

    }}))
  feather::write_feather(FDC,paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))
}else{
  FDC<-feather::read_feather(paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))
}

# http://r-sig-geo.2731867.n2.nabble.com/Extracting-coordinates-from-SpatialLinesDataFrame-object-td3212217.html
