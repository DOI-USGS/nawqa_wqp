# This file contains a first cut at partitioning whatever data are currently in
# your 1_wqpdata/out/data folder into these scenarios:
#   A (ideal) - sites with 12 samples per year, dissolved P, TP, and NO3
#   B (almost there, need more constituents) - sites with 12 samples per year, 1 or two of dissolved P, TP, NO3
#   C (almost there, need more frequency) - sites with 6 samples per year, dissolved P, TP, NO3
# (as you'll see, I've had to back the number of constituents down to 2 for A
# and C to get any sites for AL)

library(tidyverse)
library(feather) # install with install.packages('feather')
library(scipiper) # install with devtools::install_github('USGS-R/scipiper')

# pull the files off of Google Drive - browser window for authentication should
# pop up the first time you do this
indicator_files <- dir('1_wqpdata/out/data', pattern='*\\.feather.ind$', full.names=TRUE) %>%
  grep('(nitrate|phosphorus)', ., value=TRUE)
sapply(indicator_files, scipiper::gd_get, USE.NAMES=FALSE)

# read in the data and combine into one big file
data_files <- sapply(indicator_files, scipiper::as_data_file, USE.NAMES=FALSE)
combined_data <- bind_rows(lapply(data_files, function(dfile) {
  cname <- gsub('1_wqpdata/out/data/Alabama_', '', gsub("_00.\\.feather", '', dfile))
  dat <- feather::read_feather(dfile) %>%
    mutate(ConstituentGroup=cname)
}))

# need to improve the following separation into TP, dissolved P here - get
# crosswalk from Melissa. Here's a very coarse first cut:
dissolvedP <- c(
  "Inorganic phosphorus",
  "Phosphate",
  "Phosphate-phosphorus",
  "Phosphate-phosphorus as P",
  "Phosphate-phosphorus as PO4",
  "Phosphorus, phosphate (PO4) as orthophosphate",
  "Soluble Reactive Phosphorus (SRP)",
  "Polyphosphate")
totalP <- c(
  "Phosphorus, hydrolyzable plus orthophosphate",
  "Phosphorus",
  "Phosphorus as P",
  "Total Phosphorus, mixed forms")
data_split_P <- combined_data %>%
  mutate(ConstituentGroup = ifelse(
    ConstituentGroup!='phosphorus',
    ConstituentGroup,
    ifelse(CharacteristicName %in% dissolvedP, 'dissolvedP',
           ifelse(CharacteristicName %in% totalP, 'totalP',
                  NA)))) %>%
  filter(!is.na(ConstituentGroup))

# count the number of target constituents (should be TP, dissolved P, and NO3,
# but for now we'll just do any-P and NO3) per sample. this takes a few seconds,
# may get to be unweildy when we have more states and constituents, such that
# maybe breaking it up into years or states will make sense.
constituents_per_sample <- combined_data %>%
  select(ConstituentGroup, MonitoringLocationIdentifier, ResultMeasureValue, ActivityStartDateTime) %>%
  mutate(Year = as.integer(format(ActivityStartDateTime, '%Y'))) %>%
  group_by(MonitoringLocationIdentifier, ActivityStartDateTime, Year) %>%
  summarize(NumConstituentsPerSample=length(unique(ConstituentGroup))) %>%
  ungroup()

# count the number of samples per year, separating into samples with one
# constituent or two
samples_per_year <- constituents_per_sample %>%
  group_by(MonitoringLocationIdentifier, Year) %>%
  summarize(
    NumSamplesPerYear_Min1Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 1]),
    NumSamplesPerYear_Min2Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 2]),
    NumSamplesPerYear_Min3Const=length(ActivityStartDateTime[NumConstituentsPerSample >= 3])) %>%
  ungroup() %>%
  gather(temporary, NumSamplesPerYear, starts_with('NumSamplesPerYear')) %>%
  extract(temporary, 'MinNumConstituents', "NumSamplesPerYear_Min([[:alnum:]])Const") %>%
  mutate(MinNumConstituents=as.integer(MinNumConstituents))

# compute sampling frequencies for samples containing at least 1, at least 2, or
# 3 different constitutent types
sample_rates <- samples_per_year %>%
  # make sure all 5 years are represented for each site, even if there were 0
  # samples some years
  filter(Year %in% 2011:2015) %>%
  complete(Year = 2011:2015) %>%
  mutate(NumSamplesPerYear = ifelse(is.na(NumSamplesPerYear), 0, NumSamplesPerYear)) %>%
  # compute various stats for frequency over all target years
  group_by(MonitoringLocationIdentifier, MinNumConstituents) %>%
  summarize(
    MeanSamplesPerYear=mean(NumSamplesPerYear),
    MinSamplesPerYear=min(NumSamplesPerYear),
    MaxSamplesPerYear=max(NumSamplesPerYear),
    MedianSamplesPerYear=median(NumSamplesPerYear))

# split into scenarios
scenario_empty <- sample_rates %>%
  filter(MinNumConstituents == 3, MinSamplesPerYear >= 12)
scenario_A <- sample_rates %>%
  filter(MinNumConstituents == 2, MinSamplesPerYear >= 12)
scenario_B <- sample_rates %>%
  filter(MinNumConstituents == 1, MinSamplesPerYear >= 12)
scenario_C <- sample_rates %>%
  filter(MinNumConstituents == 2, MinSamplesPerYear >= 6)

#### explore results ####

# there are at least a few thousand of each, many fewer for dissolved P
data_split_P$ConstituentGroup %>% table

# site counts for each scenario - A and B are relatively infrequent, but C is
# pretty common
nrow(scenario_A)
nrow(scenario_B)
nrow(scenario_C)

# histogram of samples per year grouped by minimum number of constituents per
# sample - shows that samples with 3 constituents don't currently exist in AL,
# samples that have 1 almost always have 2 constituents, and many more sites
# have 5-10 samples per year than >= 12
sample_rates %>%
  mutate(MinNumConstituents=as.factor(MinNumConstituents)) %>%
  ggplot(aes(x=MeanSamplesPerYear, group=MinNumConstituents, fill=MinNumConstituents)) +
  geom_histogram(binwidth=1, position='dodge')

# for a map, we'll need to pull in site coordinates. this can be done with
# dataRetrieval and will eventually be done in phase 1_wqpdata
