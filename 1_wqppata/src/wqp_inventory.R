#### inventory ####

# prepare a data.frame that maps state names to the FIPS codes used by WQP
get_wqp_state_codes <- function() {
  states_xml <- xml2::read_xml('https://www.waterqualitydata.us/Codes/statecode?countrycode=US')
  states_list <- xml2::as_list(states_xml)
  states_df <- bind_rows(lapply(states_list[names(states_list)=='Code'], function(code) {
    data_frame(
      value = attr(code, 'value'),
      name = attr(code, 'desc'),
      providers = attr(code, 'providers'))
  }))
  if(nrow(states_df) < 1){stop('No State FIPS codes downloaded,
                               which is critical for later steps, check 
                               that get_wqp_state_codes works')}
  return(states_df)
}

# acquire a data frame of site x constituent information, with counts of
# observations per site-constituent combination and all the site metadata that
# looks useful
inventory_wqp <- function(ind_file, wqp_state_codes, wqp_states, wqp_codes) {
  # convert states list to FIPS list
  state_codes <- filter(wqp_state_codes, name %in% wqp_states) %>% pull(value)
  
  # identify available constituent sets
  constituents <- names(wqp_codes$characteristicName)
  
  # prepare the args to whatWQPdata. all arguments will be the same every time
  # except characteristicName, which we'll loop through to get separate counts
  # for each
  wqp_args <- list(
    statecode=state_codes,
    siteType=wqp_codes$siteType,
    characteristicName=NA, # to be updated each time through loop
    sampleMedia=wqp_codes$sampleMedia
    # we'd include dates, but they get ignored by the service behind whatWQPdata
  )
  
  # loop over the constituents, getting rows for each
  sample_time <- system.time({
    samples <- bind_rows(lapply(constituents, function(constituent) {
      message(Sys.time(), ': getting inventory for ', constituent)
      wqp_args$characteristicName <- wqp_codes$characteristicName[[constituent]]
      tryCatch({
        wqp_wdat <- do.call(whatWQPdata, wqp_args)
        mutate(wqp_wdat, constituent=constituent)
      }, error=function(e) {
        # keep going IFF the only error was that there weren't any matching sites
        if(grepl('arguments imply differing number of rows', e$message)) {
          NULL
        } else {
          stop(e)
        }
      })
    }))
  })
  message(sprintf('sample inventory complete, required %0.0f seconds', sample_time[['elapsed']]))
  
  # get additional site information
  message(Sys.time(), ': getting additional site data')
  site_time <- system.time({
    wqp_site_args <- wqp_args[names(wqp_args) != 'characteristicName']
    sites <- do.call(whatWQPsites, wqp_site_args)
  })
  message(sprintf('site inventory complete, required %0.0f seconds', site_time[['elapsed']]))
  
  # merge constituent info with site info
  wqp_info <- left_join(
    samples %>%
      select(Constituent=constituent, MonitoringLocationIdentifier, resultCount,
             MonitoringLocationName, MonitoringLocationTypeName, ResolvedMonitoringLocationTypeName),
    sites %>%
      select(MonitoringLocationIdentifier, LatitudeMeasure, LongitudeMeasure, HorizontalCoordinateReferenceSystemDatumName,
             HUCEightDigitCode, CountryCode, StateCode, CountyCode, OrganizationFormalName) %>%
      # replace lat/lon numeric flags with NAs
      mutate(LatitudeMeasure = ifelse(LatitudeMeasure < 10, NA, LatitudeMeasure),
             LongitudeMeasure = ifelse(LongitudeMeasure > -10, NA, LongitudeMeasure)),
    by='MonitoringLocationIdentifier')
  
  # write the data file and the indicator file
  data_file <- as_data_file(ind_file)
  feather::write_feather(wqp_info, path=data_file)
  gd_put(ind_file, data_file) # sc_indicate(ind_file, data_file=data_file)
  
  invisible()
}
