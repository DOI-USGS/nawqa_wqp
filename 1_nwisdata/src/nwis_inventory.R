# acquire a data frame of site x constituent information, with counts of
# observations per site-constituent combination and all the site metadata that
# looks useful
inventory_nwis <- function(
  ind_file, nwis_states=c('WI','OH','PA','NC'),
  nwis_codes=list(parameterCd=list(temperature='00010'), siteType=c('ST'), service='uv')) {

  # identify available constituent sets
  constituents <- names(nwis_codes$parameterCd)

  # prepare the args to whatWQPdata. all arguments will be the same every time
  # except characteristicName, which we'll loop through to get separate counts
  # for each
  nwis_args <- list(
    stateCd=NA, # to be updated each time through states loop below
    siteType=nwis_codes$siteType,
    service=nwis_codes$service,
    parameterCd=NA # to be updated each time through constituents loop below
    # we'd include startDt and endDt, but they get ignored by the service behind whatNWISdata
  )

  # loop over the constituents and states, getting rows for each
  sample_time <- system.time({
    samples <- bind_rows(lapply(constituents, function(constituent) {
      nwis_args$parameterCd <- nwis_codes$parameterCd[[constituent]]
      samples_1state <- bind_rows(lapply(nwis_states, function(state) {
        message(sprintf('%s: getting inventory for %s in %s', Sys.time(), constituent, state))
        nwis_args$stateCd <- state
        tryCatch({
          nwis_wdat <- do.call(whatNWISdata, nwis_args)
          mutate(
            nwis_wdat,
            constituent=constituent,
            dec_lat_va=as.numeric(dec_lat_va),
            dec_long_va=as.numeric(dec_long_va),
            alt_acy_va=as.character(alt_acy_va),
            begin_date=as.Date(begin_date),
            end_date=as.Date(end_date))
        }, error=function(e) {
          # keep going IFF the only error was that there weren't any matching sites
          # if(grepl('arguments imply differing number of rows', e$message)) {
          #   NULL
          # } else {
            stop(e)
          # }
        })
      }))
    }))
  })
  message(sprintf('sample inventory complete, required %0.0f seconds', sample_time[['elapsed']]))

  ## NOTE: the count_nu field returned by whatNWISdata is simply end_date - ##
  ## begin_date. Confirmed by Brad Garner on 7/26/2018 ##

  # we could call whatNWISsites next, but all the columns provided by that query
  # are also provided by the above inventory query already

  # write the data file and the indicator file. use RDS to retain the attributes
  data_file <- as_data_file(ind_file)
  saveRDS(samples, file=data_file)
  gd_put(ind_file, data_file) # sc_indicate(ind_file, data_file=data_file)

  invisible()
}
