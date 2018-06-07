#### pull ####

# we've been asked to please NOT use a cluster because a postdoc took down the
# entire WQP a couple of years ago by doing that, and 12 simultaneous requests
# would be way too many. Alternative ways we can speed up the pull:
# * subset spatially rather than temporally - loop over state/HUC rather than years
# * probably fastest to pull all variables at once
# * probably faster to request whole states at once rather than giving explicit site lists

# wrap dplyr::filter for use within a remake command
filter_partitions <- function(partitions, pull_task) {
  dplyr::filter(partitions, PullTask==pull_task)
}

# package all the configuration information for a single partition into each row
# so that we know whether to re-pull a partition based solely on whether the
# contents of that row have changed. each row should therefore include the lists
# of sites and characteristic names.
partition_inventory <- function(ind_file, inventory_ind, wqp_pull, wqp_state_codes, wqp_codes) {
  # read in the inventory, which includes all dates and is therefore a superset
  # of what we'll be pulling for a specific date range

  #scipiper::sc_retrieve(inventory_ind)
  inventory <- feather::read_feather(scipiper::as_data_file(inventory_ind)) %>%
    select(Constituent, Site=MonitoringLocationIdentifier, StateCode, resultCount)

  # split up the work into manageably sized partitions (~max size specified in wqp_pull)
  partitions <- bind_rows(lapply(unique(inventory$Constituent), function(constituent) {
    bind_rows(lapply(unique(inventory$StateCode), function(statecode) {
      # an atomic group is a combination of parameters that can't be reasonably
      # split into multiple WQP pulls - in this case we're defining atomic
      # groups as distinct combinations of constituent (a group of
      # characteristicNames) and site ID. We could potentially split
      # constituents further, or pull multiple constituents at once, or pull an
      # entire state at once, but we're doing it this way
      atomic_groups <- inventory %>%
        filter(StateCode==statecode, Constituent == constituent) %>%
        group_by(Site) %>%
        summarize(NumObs=sum(resultCount)) %>%
        arrange(desc(NumObs))

      # split the full pull (combine atomic groups) into right-sized partitions
      # by an old but fairly effective paritioning heuristic: pick the number of
      # partitions desired, sort the atomic groups by descending size, and then
      # go down the list, each time adding the next atomic group to the
      # partition that's currently smallest
      num_partitions <- ceiling(sum(atomic_groups$NumObs) / wqp_pull$target_pull_size)
      partition_sizes <- rep(0, num_partitions)
      assignments <- rep(0, nrow(atomic_groups))
      for(i in seq_len(nrow(atomic_groups))) {
        size_i <- atomic_groups[[i,"NumObs"]]
        smallest_partition <- which.min(partition_sizes)
        assignments[i] <- smallest_partition
        partition_sizes[smallest_partition] <- partition_sizes[smallest_partition] + size_i
      }

      # create and return data.frame rows for binding in the two lapply loops
      state <- wqp_state_codes %>%
        filter(value==sprintf('US:%s', statecode)) %>%
        pull(name)
      atomic_groups %>%
        mutate(
          State=state,
          StateCode=statecode,
          Constituent=constituent,
          NumPulls=num_partitions,
          PullTask=sprintf('%s_%s_%03d', gsub(' ', '_', state), constituent, assignments))
    }))
  }))

  # prepare nested collections of other possible arguments to include in the readWQPdata calls
  parameter_codes <- bind_rows(lapply(unique(partitions$Constituent), function(constituent) {
    data_frame(Constituent=constituent, CharacteristicNames=wqp_codes$characteristicName[[constituent]]) %>%
      tidyr::nest(-Constituent, .key='Params')
  }))
  site_types <- data_frame(SiteType=wqp_codes$siteType)
  pull_tasks_df <- partitions %>%
    dplyr::group_by(PullTask, NumPulls, State, StateCode, Constituent) %>%
    tidyr::nest(Site, NumObs, .key="Sites") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(parameter_codes, by='Constituent') %>%
    dplyr::mutate(SiteType=list(site_types))

  # write the data file and the indicator file
  data_file <- as_data_file(ind_file)
  saveRDS(pull_tasks_df, file=data_file)
  gd_put(ind_file, data_file) # sc_indicate(ind_file, data_file=data_file)

  invisible()
}

plan_wqp_pull_per_constituent <- function(constituents, folders, folders_item) {

  partitions <- scipiper::create_task_step(
    step_name = 'partitions',
    target_name = function(task_name, step_name, ...) {
      sprintf('partitions_%s', task_name)
    },
    command = function(task_name, ...) {
      sprintf("readRDS('%s')",
              file.path(folders$tmp, sprintf("partition_%s.rds", task_name)))
    }
  )

  make_plan <- scipiper::create_task_step(
    step_name = 'make_plan',
    target_name = function(task_name, step_name, ...) {
      sprintf('wqp_pull_plan_%s', task_name)
    },
    command = function(task_name, ...) {
      sprintf("plan_wqp_pull(%s, I('%s'), %s)",
              sprintf('partitions_%s', task_name),
              task_name,
              folders_item)
    }
  )

  create_plan <-  scipiper::create_task_step(
    step_name = 'create_plan',
    target_name = function(task_name, step_name, ...) {
      sprintf('tasks_1_wqp_%s.yml', task_name)
    },
    command = function(task_name, ...) {
      sprintf(
        paste(
          "create_task_makefile(",
          "makefile=target_name,",
          "task_plan=%s,",
          "include=I(c('1_wqpdata.yml',target_name)),",
          "packages=I(c('dplyr', 'dataRetrieval', 'feather', 'scipiper')),",
          "file_extensions=I(c('ind','feather')))",
          sep="\n      "
        ),
        sprintf('wqp_pull_plan_%s', task_name)
      )
    }
  )

  pull_data <-  scipiper::create_task_step(
    step_name = 'pull_data',
    target_name = function(task_name, step_name, ...) {
      sprintf('1_wqpdata/log/tasks_1_wqp_%s.ind', task_name)
    },
    command = function(task_name, ...) {
      sprintf(
        paste(
          "loop_tasks(",
          "task_plan=%s,",
          "task_makefile='%s',",
          "num_tries=I(30), sleep_on_error=I(20))",
          sep="\n      "
        ),
        sprintf('wqp_pull_plan_%s', task_name),
        sprintf('tasks_1_wqp_%s.yml', task_name)
      )
    }
  )

  task_plan <- scipiper::create_task_plan(
    task_names=sort(names(constituents$characteristicName)),
    task_steps=list(partitions, make_plan, create_plan, pull_data),
    final_steps='pull_data',
    add_complete=FALSE,
    ind_dir=folders$log)

}

# prepare a plan for downloading (from WQP) and posting (to GD) one data file
# per state per constituent
plan_wqp_pull <- function(partitions, constituent, folders) {

  partition <- scipiper::create_task_step(
    step_name = 'partition',
    target_name = function(task_name, step_name, ...) {
      sprintf('partition_%s', task_name)
    },
    command = function(task_name, ...) {
      sprintf("filter_partitions(partitions_%s, I('%s'))", constituent, task_name)
    }
  )

  # steps: download, post
  download <- scipiper::create_task_step(
    step_name = 'download',
    target_name = function(task_name, step_name, ...) {
      scipiper::as_ind_file(file.path(folders$tmp, sprintf('%s.feather', task_name)))
    },
    command = function(task_name, ...) {
      paste(
        "get_wqp_data(",
        "ind_file=target_name,",
        sprintf("partition=partition_%s,", task_name),
        "wq_dates=wq_dates)",
        sep="\n      ")
    }
  )

  post <- scipiper::create_task_step(
    step_name = 'post',
    target_name = function(task_name, step_name, ...) {
      scipiper::as_ind_file(file.path(folders$out, sprintf('%s.feather', task_name)))
    },
    command = function(task_name, ...) {
      sprintf(
        paste(
          "gd_put(",
          "remote_ind=target_name,",
          "local_source='%s',",
          "mock_get=I('move'),",
          "on_exists=I('update'))",
          sep="\n      "),
        scipiper::as_ind_file(file.path(folders$tmp, sprintf('%s.feather', task_name))))
    }
  )

  retrieve <- scipiper::create_task_step(
    step_name = 'retrieve',
    target_name = function(task_name, step_name, ...) {
      file.path(folders$out, sprintf('%s.feather', task_name))
    },
    command = function(task_name, target_name, ...) {
      sprintf(
        paste(
          "gd_get(",
          "ind_file='%s')",
          sep="\n      "),
        scipiper::as_ind_file(target_name))
    }
  )

  task_plan <- scipiper::create_task_plan(
    task_names=sort(partitions$PullTask),
    task_steps=list(partition, download, post, retrieve),
    final_steps='post',
    add_complete=FALSE,
    ind_dir=folders$log)

}

get_wqp_data <- function(ind_file, partition, wq_dates) {

  # prepare the argument to pass to readWQPdata
  wqp_args <- list(
    characteristicName=partition$Params[[1]]$CharacteristicNames,
    startDateLo=wq_dates$startDate,
    startDateHi=wq_dates$endDate
  )
  # specify sites either by whole state and siteType (if we didn't split the
  # site into multiple site sets) or as a specific list of sites (otherwise)
  if(partition$NumPulls == 1) {
    wqp_args <- c(wqp_args, list(
      statecode=partition$StateCode,
      siteType=partition$SiteType[[1]]$SiteType
    ))
  } else {
    wqp_args <- c(wqp_args, list(
      siteid=partition$Sites[[1]]$Site
    ))
  }

  # do the data pull
  wqp_dat_time <- system.time({
    wqp_dat <- do.call(dataRetrieval::readWQPdata, wqp_args)
  })
  message(
    "WQP pull for ", partition$PullTask,
    " took ", wqp_dat_time['elapsed'], " seconds and",
    " returned ", nrow(wqp_dat), " rows")

  # make wqp_dat a tibble, converting either from data.frame (the usual case) or
  # NULL (if there are no results)
  wqp_dat <- as_data_frame(wqp_dat)

  # write the data and indicator file. do this even if there were 0 results
  # because remake expects this function to always create the target file
  data_file <- as_data_file(ind_file)
  feather::write_feather(as_data_frame(wqp_dat), path=data_file)
  sc_indicate(ind_file, data_file=data_file)

  invisible()
}
