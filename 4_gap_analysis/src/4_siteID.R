# Download data for streamgages currently active that we can couple with WQ Monitoring sites
# Warning - this will take a long time to run

library(tidyverse)
library(feather)
library(XML)
library(xml2)

# Read in URL from NWIS Service
xml.url<-readChar("4_gap_analysis/in/XML.txt",file.info("4_gap_analysis/in/XML.txt")$size)

# Read in XML
xmlLines<-read_xml(xml.url)
# First hierarchy level
xmlObs<-xmlLines %>% xml_find_all("site")
# Find site number
STAID<-xmlObs %>% xml_find_all("site_no") %>% xml_text() %>% as.character()
# Get Lat and Long
LAT<-xmlObs %>% xml_find_all("dec_lat_va") %>% xml_text() %>% as.character()
LONG<-xmlObs %>% xml_find_all("dec_long_va") %>% xml_text() %>% as.character()

# for each station, get entire POR, write out data, and get number of data
countQ<-function(station){
  Q<-dataRetrieval::readNWISdata(site=station,service="dv",asDateTime=TRUE,parameterCd="00060",
                                 startDate="1900-01-01",endDate="2015-12-31")
  write.csv(Q,paste("Q/",station,".csv",sep=""),quote=F,row.names=F)
  numDays<-length(!is.na(Q$X_00060_00003))
  return(numDays)
}

QRec<-lapply(STAID,countQ)

# File of gages to use
gages<-data.frame(STA=STAID,X=LONG,Y=LAT,days=unlist(QRec))
write_feather(gages,"gages.feather")

