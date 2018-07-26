count_wtemp_sites <- function(nwis_inv_ind, wqp_inv_ind) {
  # Use the NWIS inventory to learn about continuous water temperature data.
  # Eventually it'd be nice to pull down all the NWIS samples, but until then,
  # we need to take shortcuts to approximately determine which sites have good
  # continuous monitoring coverage
  nwis_inventory <- readRDS(sc_retrieve(nwis_inv_ind))
  continuous <- nwis_inventory %>%
    filter(end_date - begin_date >= as.difftime(3000, units='days')) %>% # arbitrary choice to define "continuous" as sites having monitoring occurring over a 3000-day range
    pull(site_no)

  # Temporarily use the WQP inventory to learn about historical water
  # temperature data, but we'll want to switch to counting actual samples and
  # their dates to set good criteria for whether a site has historical coverage
  # or not
  wqp_inventory <- feather::read_feather(sc_retrieve(wqp_inv_ind))
  historical <- wqp_inventory %>%
    filter(Constituent == 'temperature') %>%
    filter(resultCount > 500) %>% # placeholder choice (needs revision) for how to decide whether a site has historical coverage
    mutate(site_no = gsub('USGS-', '', MonitoringLocationIdentifier)) %>%
    pull(site_no)

  # Count and return the number of sites that have both continuous and
  # historical coverage
  continuous_w_history <- intersect(continuous, historical) %>% length()
  return(continuous_w_history)
}
