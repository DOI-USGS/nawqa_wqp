library(dataRetrieval)

sapply(list.files("1_wqpdata/src",full.names = T),source,.GlobalEnv)

staID<-unique(constituents_per_sample$MonitoringLocationIdentifier)

test<-dataRetrieval::readWQPdata(staID[1],"Phosphorus")
