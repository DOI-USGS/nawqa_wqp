# This identifies closest streamgage and NHDflowline to each WQM site

library(tidyverse)
library(rgeos)
library(RANN)

# Select the constituent Group
Con<-"Nut"

# Read in files
# NHD Directory, lists the regional directorys
nhdDir<-list.dirs("4_gap_analysis/in/GIS",recursive=F)
# Site Files by state
siteFiles<-dir(paste("4_gap_analysis/in/Site/",Con,sep=""), pattern="*sites.feather",full.names = T)
# Gages derived from an XML file
gages<-feather::read_feather("gages.feather")

# Open site files and bind together in a data frame for a given Con or Con group
siteInfo <- bind_rows(
  lapply(siteFiles, function(dfile) {
    dat <- feather::read_feather(dfile)
  })
)

# Eliminate duplicate sites
siteInfo<-siteInfo %>% distinct(MonitoringLocationIdentifier,.keep_all=T)
# Set row names as site IDs
rownames(siteInfo)<-siteInfo$MonitoringLocationIdentifier
# Keep only the lat and long columns
siteInfo2<-siteInfo %>% select(dec_lon_va,dec_lat_va)
colnames(siteInfo2) <- c("lon","lat")

# convert to albers for proximity analysis and distance units (meters)
sp::coordinates(siteInfo2)<- ~lon+lat
# define projection as WGS84
sp::proj4string(siteInfo2)<-sp::CRS("+init=epsg:4269")
# Re-project to albers
siteInfo_alb <- sp::spTransform(siteInfo2, sp::CRS("+init=epsg:5070"))

# Convert gages to albers with same procedure for the sites
gages$X<-as.numeric(as.character(gages$X))
gages$Y<-as.numeric(as.character(gages$Y))
gages$STA<-as.character(gages$STA)
colnames(gages)<-c("STA","lon","lat","days")
gages<-gages[!is.na(gages$lon),]
gagesXY <- gages %>% select(lon,lat)
rownames(gagesXY)<-gages$STA
sp::coordinates(gagesXY)<- ~lon+lat
sp::proj4string(gagesXY)<-sp::CRS("+init=epsg:4269")
# Re-project to albers
gagesAlb <- sp::spTransform(gagesXY, sp::CRS("+init=epsg:5070"))
rownames(gagesAlb@coords)<-gages$STA

# Run nearest neighbor analysis on the gages and sample sites
closestSTAID<-as.data.frame(nn2(gagesAlb@coords,siteInfo_alb@coords,k=3))

# Bind 3 closest staions and their distances (meters) into a data frame
STAIDs<-as.data.frame(cbind(rownames(gagesAlb@coords)[closestSTAID[,1]],
                            rownames(gagesAlb@coords)[closestSTAID[,2]],
                            rownames(gagesAlb@coords)[closestSTAID[,3]],
                            closestSTAID[,4],closestSTAID[,5],closestSTAID[,6]))
STAIDs$monID<-rownames(siteInfo)
rownames(STAIDs)<-rownames(siteInfo)
feather::write_feather(STAIDs,paste("4_gap_analysis/in/Site/",Con,"/gagesdist_Alb.feather",sep=""))

# Run the same exercise but with NHDPlus flowlines
if (!file.exists(paste0("4_gap_analysis/in/Site/",Con,"/NHDPlus01_alb_index.feather"))){
  for (nhd in nhdDir){
    print(nhd)
    flowLines<-rgdal::readOGR(paste(nhd[1],"/NHDSnapshot/Hydrography",sep=""),"NHDFlowline")
    # Re-project to albers
    flowLines_alb <- sp::spTransform(flowLines, sp::CRS("+init=epsg:5070"))
    # COMID can also be ComID depending on the region (just changing by hand for now)
    IDs<-flowLines_alb$COMID

    # XY points for each line
    res <- lapply(slot(flowLines_alb, "lines"), function(x) lapply(slot(x, "Lines"),
                                                               function(y) slot(y, "coords")))
    # Flatten point XYs into a list
    res2<-flatten(res)
    names(res2)<-IDs
    dFrame<-do.call(rbind,res2)
    # This is getting the three closest points, not three closest segments
    res_names<-Map(cbind,res2,ID=names(res2))
    dFrame_names<-do.call(rbind,res_names)
    closest<-nn2(dFrame,siteInfo_alb@coords,k=1)

    comIDs<-as.data.frame(cbind(dFrame_names[closest$nn.idx[,1],3]))

    # Write out the COMIDs
    fileName<-unlist(strsplit(nhd[1],"/"))[5]
    feather::write_feather(comIDs,paste("4_gap_analysis/in/Site/",Con,"/",fileName,"_alb_index.feather",sep=""))
    # write out the distances
    distDF<-as.data.frame(closest$nn.dists)
    feather::write_feather(distDF,paste("4_gap_analysis/in/Site/",Con,"/",fileName,"_alb_value.feather",sep=""))
  }
}

# Open site files
nhdFiles<-list.files(paste0("4_gap_analysis/in/Site/",Con),"*_alb_index.feather",full.names=T)
distFiles<-list.files(paste0("4_gap_analysis/in/Site/",Con),"_alb_value.feather",full.names=T)
# Load NHDPlus characteristics to get the drainage area
NHDPlusChar<-read.csv("4_gap_analysis/in/BASIN_CHAR_TOT_CONUS.txt")

# Matrix of distance from site to NHDPlus feature in each region
distInfo <- data.frame(bind_cols(lapply(distFiles, function(dfile) {dat <- feather::read_feather(dfile)})))
# Get the index of the minimum feature
distIndex<-apply(distInfo,1,which.min)
# Index for which nhdplus region is closest
nhdInfo <- data.frame(bind_cols(lapply(nhdFiles, function(dfile) {dat <- feather::read_feather(dfile)})))

# Build information into a dataframe, Join with NHDPlus data to get drainage attributes
nhdClose<-data.frame(comID=as.integer(nhdInfo[cbind(1:nrow(nhdInfo), distIndex)]),region=nhdDir[c(distIndex)],
                     dist=distInfo[cbind(1:nrow(distInfo),distIndex)]) %>%
  left_join(NHDPlusChar,by=c("comID"="COMID")) %>% 'rownames<-'(rownames(siteInfo)) %>% mutate(monID=rownames(siteInfo)) %>%
  feather::write_feather(paste0("4_gap_analysis/in/Site/",Con,"/Nat_ComIDS"))

# Join the NHDPlus and gage proximity datasets
staClose<-nhdClose %>% left_join(STAIDs,by="monID") %>% select (c(comID,region,dist,TOT_BASIN_AREA,monID,V1,V2,V3,V4,V5,V6)) %>%
  mutate_if(is.factor, as.character)

# Retrieve the drainage areas for each gaging station, and convert to square kilometers
sta_DA<-dataRetrieval::readNWISsite(gages$STA) %>% mutate(DA_KM= drain_area_va * 2.589) %>% select(site_no,DA_KM)

# Difference straight magnitude or percent
sta_DA$Dist<-NA

# retrieve gage with the closest DA to flowlines in common
getSTA<-function(i){
  print (i)
  # Get the closest gages
  gagesOne<-as.vector(staClose[i,c(6:8)])
  # Get the drainage area of the COMID
  DA_comID<-staClose[i,4]
  gages_DA<-sta_DA[which(sta_DA$site_no %in% gagesOne),]
  gages_DA$Dist<-abs(gages_DA$DA_KM-DA_comID)
  gages_DA$comID<-staClose[i,]$comID
  #daysPOR<-gages[gages$STA %in% unlist(gagesOne),]
  #staList<-append(staList,gages_DA$site_no[which.min(abs(gages_DA$DA_KM-DA_comID))])
  return(gages_DA[which.min(gages_DA$Dist),])
}
gages<-bind_rows(lapply(c(1:nrow(staClose)),getSTA))
feather::write_feather(gages,paste("4_gap_analysis/in/Site/",Con,"/ClosestGage.feather",sep=""))






# if (!file.exists(paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))){
#   options(warn=1)
#   # combined data may not work
#   FDC <- bind_rows(
#     lapply(1:length(gages$STA), function(i) {
#       gage<-gages$STA[i]
#       #print(i)
#       if (file.exists(paste("Q/",gage,".csv",sep=""))){
#         # read in Q file
#       inQ<-read.csv(paste("Q/",gage,".csv",sep=""),colClasses="character")
#         # Test for non-regular headings, convert if necessary
#           if("X_.Final._00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_.Final._00060_00003}
#           if("X_.Regression._00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_.Regression._00060_00003}
#           if("X_.Primary.Stream.Flow._00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_.Primary.Stream.Flow._00060_00003}
#           if("X_Red.River.Only_00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_Red.River.Only_00060_00003}
#           if("X_FROM.DCP_00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_FROM.DCP_00060_00003}
#           if("X_Powerhouse.Releases_00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_Powerhouse.Releases_00060_00003}
#           if("X_.Data.from.10.1.1992.Forward._00060_00003" %in% colnames(inQ)){
#             inQ$X_00060_00003<-
#               inQ$X_.Data.from.10.1.1992.Forward._00060_00003}
#
#         inQ<-inQ[c("agency_cd","site_no","dateTime","X_00060_00003")]
#         inQ$dateTime<-as.Date(inQ$dateTime)
#         inQ$X_00060_00003<-as.numeric(inQ$X_00060_00003)
#         inQ$WaterYear<-ifelse(month(inQ$dateTime)<10,year(inQ$dateTime),year(inQ$dateTime)+1)
#         # get Q from each station for each date
#         # 3 gages have NAs for entire time series and will generate soft warnings
#         # stats
#         # 1 - FDC (entire POR); 2-FDC (by year); 3 - 7-day MA; 4 - 15-day MA; 5 - 31-day MA;
#         # 6 - 7-day rank; 7 - 15-day rank; 8 - 31-day rank
#         nearQ_singleGage<-inQ %>% filter(site_no==as.character(gage)) %>% mutate(fdr=hydroTSM::fdc(X_00060_00003,plot=F)) %>%
#           group_by(WaterYear) %>% mutate(fdrAnn=hydroTSM::fdc(X_00060_00003,plot=F)) %>% ungroup() %>%
#           mutate(ma7=rollmean(X_00060_00003,7,malign='center',fill=NA),ma15=rollmean(X_00060_00003,15,malign='center',fill=NA),
#                  ma31=rollmean(X_00060_00003,31,malign='center',fill=NA)) %>%
#           mutate(rank7=zoo::rollapply(inQ$X_00060_00003,15,rank,fill=NA,align="center"),
#                  rank15=zoo::rollapply(inQ$X_00060_00003,15,rank,fill=NA,align="center"),
#                  rank31=zoo::rollapply(inQ$X_00060_00003,31,rank,fill=NA,align="center"))
#     }}))
#   feather::write_feather(FDC,paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))
# }else{
#   FDC<-feather::read_feather(paste("4_gap_analysis/in/Q/",Con,"/FDR.feather",sep=""))
# }

# http://r-sig-geo.2731867.n2.nabble.com/Extracting-coordinates-from-SpatialLinesDataFrame-object-td3212217.html
