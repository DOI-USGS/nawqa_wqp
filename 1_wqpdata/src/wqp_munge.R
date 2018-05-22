# this file requires DoMC and foreach packges

#### munge & merge tasks ####

plan_wqp_munge <- function(partitions, pull_plan, folders) {

  # We will create 1 file per constituent (all states) - should be < 1.5 Gb per
  # file for these constituents. Might want to create more, smaller files if
  # adding constituents with more observations
  constituents <- unique(partitions$Constituent)

  # Extract the raw-file target names of the WQP pull plan into a vector named by
  # the PullTask, which is also given in the partitions table
  raw_files <- sapply(pull_plan, function(task) { task$steps$retrieve$target_name })

  combine <- scipiper::create_task_step(
    step_name = 'combine',
    target_name = function(task_name, ...) {
      file.path(folders$tmp, sprintf('all_raw_%s.feather', task_name))
    },
    command = function(task_name, ...) {
      paste(
        "combine_feathers(",
        sprintf("data_file=target_name,"),
        partitions %>%
          filter(Constituent==task_name) %>%
          pull(PullTask) %>%
          raw_files[.] %>%
          sprintf("'%s'", .) %>%
          paste(collapse=",\n      ") %>%
          paste0(")"),
        sep="\n      ")
    }
  )

  munge <- scipiper::create_task_step(
    step_name = 'munge',
    target_name = function(task_name, ...) {
      file.path(folders$tmp, sprintf('all_%s.feather', task_name))
    },
    command = function(task_name, steps, ...) {
      paste(
        sprintf("munge_%s(", task_name),
        "data_file=target_name,",
        sprintf("raw_file='%s')", steps$combine$target_name),
        sep="\n      ")
    }
  )

  post <- scipiper::create_task_step(
    step_name = 'post',
    target_name = function(task_name, ...) {
      scipiper::as_ind_file(file.path(folders$out, sprintf('all_%s.feather', task_name)))
    },
    command = function(steps, ...) {
      paste(
        "gd_put(",
        "remote_ind=target_name,",
        sprintf("local_source='%s',", steps$munge$target_name),
        "mock_get=I('copy'),",
        "on_exists=I('update'))",
        sep="\n      ")
    }
  )

  task_plan <- scipiper::create_task_plan(
    task_names=sort(constituents),
    task_steps=list(combine, munge, post),
    final_steps='post',
    add_complete=FALSE,
    ind_dir=folders$log)
}

create_wqp_munge_makefile <- function(makefile, task_plan, pull_makefile) {
  create_task_makefile(
    makefile=makefile, task_plan=task_plan,
    include=pull_makefile, # include the pull_makefile (which includes remake.yml) so we know how to build the raw data files if needed
    packages=c('tidyverse', 'feather', 'scipiper','foreach','doMC'),
    file_extensions=c('ind','feather'))
}

#### merge ####

combine_feathers <- function(data_file, ...) {
  # read and combine the individual raw data pull files
  feathers <- c(...)

  # most columns are not particularly useful and they can be different by state so select
  #only the vital ones here:
  unified.cols <- c(
    'ActivityStartDate','ResultMeasureValue','ResultMeasure.MeasureUnitCode'
    ,'MonitoringLocationIdentifier','CharacteristicName','OrganizationFormalName',
    'OrganizationIdentifier','ActivityStartTime.Time')

  #Use foreach and DoMC to make this operation parallel for computers
  #with more than 3 cores

  cores <- detectCores()-2
  if(cores < 1){cores <- 1}
  registerDoMC(cores=cores)
  combo <- foreach(i=1:length(feathers),.combine=rbind) %dopar% {
    library(feather)
    feather::read_feather(feathers[i],columns=unified.cols)
  }


  # write the data file
  feather::write_feather(combo, path=data_file)
}

#### munge ####

# munge raw_file to data_file for each constituent
munge_doc <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_doc.feather') {
  raw_df <- feather::read_feather(raw_file)
  # raw_df %>% group_by(ResultMeasure.MeasureUnitCode) %>%
  #   summarize(n=n()) %>%
  #   View(.)

  #DOC data is quite messy, might only keep mg/L for now. Lots of % which is confusing
  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}


munge_poc <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_poc.feather') {
  raw_df <- feather::read_feather(raw_file)
  # raw_df %>% group_by(ResultMeasure.MeasureUnitCode) %>%
  #   summarize(n=n()) %>%
  #   View(.)

  #POC data is quite messy, might only keep mg/L for now. Lots of % which is confusing
  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}


munge_cdom <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_cdom.feather') {
  raw_df <- feather::read_feather(raw_file)

  # > t(t(table(gsub(' ', '', raw_df$ResultMeasure.MeasureUnitCode), useNA='always')))
  # mg/l        3
  # None       33
  # RFU       673
  # ug/l QSE  568
  # <NA>        4

  # preserving the original units, at least for now
  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}

munge_nitrate <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_nitrate.feather') {
  raw_df <- feather::read_feather(raw_file)


  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}

munge_tn <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_tn.feather') {
  raw_df <- feather::read_feather(raw_file)
  # raw_df %>% group_by(ResultMeasure.MeasureUnitCode) %>%
  #   summarize(n=n()) %>%
  #   View(.)

  #DOC data is quite messy, might only keep mg/L for now. Lots of % which is confusing
  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}

munge_p <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_p.feather') {
  raw_df <- feather::read_feather(raw_file)


  unit_map <- NA

  munge_any(raw_df, data_file, unit_map)
}


munge_chlorophyll <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_chlorophyll.feather') {
  raw_df <- feather::read_feather(raw_file)

  # > t(t(table(gsub(' ', '', raw_df$ResultMeasure.MeasureUnitCode), useNA='always')))
  #           28704
  # #/100ml       1
  # %           504
  # g/m2         42
  # IVFU         90
  # mg           52
  # mg/cm3       67
  # mg/l       9880
  # mg/m2     36999
  # mg/m3    369130
  # ml           15
  # None      29603
  # NTU           5
  # ppb       40008
  # ppm         805
  # RFU        3207
  # ug/cm2      244
  # ug/l    2064753
  # volts      5083
  # <NA>     161634

  # convert volumetric to ug/l, leave everything else alone. this may or may not
  # be the best approach - matt, please weigh in
  unit_map <- data_frame(
    UnitsRaw=c('mg/l', 'mg/cm3', 'mg/m3', 'ppb', 'ppm', 'ug/l'),
    Conversion = c(1000, 1000000, 1, 1, 1000, 1), # ppb and ppm conversions are approximate
    Units='ug/l')

  munge_any(raw_df, data_file, unit_map)
}

munge_secchi <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_secchi.feather'){
  raw_df <- feather::read_feather(raw_file)

  # > t(t(table(gsub(' ', '', raw_df$ResultMeasure.MeasureUnitCode), useNA='always')))
  #                361
  # cm       19903
  # degC         1
  # degF         7
  # ft      357901
  # ft/sec      20
  # in       49481
  # m      1366724
  # mg         171
  # mi           1
  # None       642
  # <NA>     42673

  # convert values with clear units to meters. 'degC', 'degF', 'ft/sec', 'mg'
  # would be candidates for exclusion, but let's leave that to @matt for final
  # decision
  unit_map <- data_frame(
    UnitsRaw=c('m', 'in', 'ft', 'cm', NA),
    Conversion = c(1, 0.0254, 0.3048, 0.01, NA),
    Units=c(rep('m', 4), NA))

  munge_any(raw_df, data_file, unit_map)
}

munge_tss <- function(data_file, raw_file='1_wqpdata/tmp/wqp/all_raw_tss.feather'){
  raw_df <- feather::read_feather(raw_file)

  # > t(t(table(gsub(' ', '', raw_df$ResultMeasure.MeasureUnitCode), useNA='always')))
  #               33
  # %           4523
  # kg            29
  # mg/l     2475410
  # None          16
  # NTU          272
  # ppm          644
  # tons/day     529
  # ug/l         456
  # <NA>      223224

  # convert values with clear units to meters, and convert values with
  # unclear/incorrect units to NA
  unit_map <- data_frame(
    UnitsRaw=c('mg/l', 'ppm', 'ug/l', NA),
    Conversion = c(1, 1, 0.001, NA),
    Units=c(rep('mg/l', 3), NA))

  munge_any(raw_df, data_file, unit_map)
}

# raw_df is the read-in data file. data_file is the filename to write to.
# unit_map should have columns UnitsRaw, Conversion, and Units (where Units is
# the final units after multiplying by Conversion). If unit_map is NA, raw units
# will be preserved.
munge_any <- function(raw_df, data_file, unit_map=NA) {
  stopifnot(is.na(unit_map) || isTRUE(is.data.frame(unit_map) && all.equal(names(unit_map), c('UnitsRaw','Conversion','Units'))))

  munged_df <- raw_df %>%
    select(
      Date=ActivityStartDate,
      Time=ActivityStartTime.Time,
      Value=ResultMeasureValue,
      UnitsRaw=ResultMeasure.MeasureUnitCode,
      SiteID=MonitoringLocationIdentifier,
      CharacteristicName=CharacteristicName # there are several methods fields, but they're all super noisy
    ) %>%
    mutate(UnitsRaw = gsub(' ', '', UnitsRaw))

  if(is.data.frame(unit_map)) {
    # implement conversions from unit_map, preserving any unknown units
    all_units <-
      data_frame(
        UnitsRaw = setdiff(unique(munged_df$UnitsRaw), unit_map$UnitsRaw),
        Conversion = 1,
        Units = UnitsRaw) %>%
      bind_rows(unit_map)
    munged_df <- munged_df %>%
      left_join(all_units, by='UnitsRaw') %>%
      mutate(Value = Value * Conversion)
  } else {
    munged_df <- munged_df %>%
      mutate(Units = UnitsRaw)
  }
  # t(t(table(gsub(' ', '', munged_df$Units), useNA='always')))

  munged_df <- munged_df %>%
    filter(!is.na(Value)) %>%
    select(Date, SiteID, Value, Units, CharacteristicName)

  feather::write_feather(munged_df, data_file)
}
