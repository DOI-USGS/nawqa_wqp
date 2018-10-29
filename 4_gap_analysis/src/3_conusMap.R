# Bugs in tidyverse interfere with ggplot2 maps functionality
# https://stackoverflow.com/questions/45066628/cannot-run-map-datastate
# library(tidyverse)
library(leaflet)
library(tidyverse)
library(mapview)
library(htmlwidgets)
library(webshot)

Con<-"WQ"

siteFiles<-dir(paste('4_gap_analysis/in/Site/',Con,sep=""), pattern='*sites.feather',full.names = T)

# combine data
siteInfo <- bind_rows(
  lapply(siteFiles, function(dfile) {
    dat <- feather::read_feather(dfile)
  })
)

siteInfo<-siteInfo %>% distinct(MonitoringLocationIdentifier,.keep_all=T) %>% filter(dec_lon_va != 0 | dec_lat_va !=0) %>%
  filter(!is.na(Tier))
# Extract character and integer from Tier
siteInfo$TierNum <- as.numeric(str_extract(siteInfo$Tier, "[0-9]+"))
siteInfo$Char<-str_extract(siteInfo$Tier, "[aA-zZ]+")
siteInfo$TierNum<-factor(siteInfo$TierNum)

# Make a simple map
m=leaflet()
m=addTiles(m)
colorFactors = colorFactor(c('blue','green','red'), # Change for tiers
                           domain = siteInfo$TierNum) #,
m = addCircleMarkers(m,
                     lng = siteInfo$dec_lon_va, # we feed the longitude coordinates
                     lat = siteInfo$dec_lat_va,
                     popup = siteInfo$Tier,#OrganizationIdentifier,
                     radius = 3,
                     stroke = FALSE,
                     fillOpacity = 0.75,
                     col=colorFactors(siteInfo$TierNum),
                     group = "1 - graffiti & noise"
) %>% addLegend(title = "Nitrate, dissolved phosphorus, and total phosphorus <br> Monthly sampling at WQP sites, 2010-2015",
                  "bottomright",pal=colorFactors,values=siteInfo$TierNum) %>%
  setView(-96.5795,39.8283,zoom=3.3)

#saveWidget(m, "d:/abock/temp/NWQ_seasonsmap.html", selfcontained = TRUE)



