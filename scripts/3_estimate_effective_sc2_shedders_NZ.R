# This script takes the data shared by Hewitt et al from the SM of their pub
# https://doi.org/10.1016/j.watres.2021.118032 and onset dates that they
# kindly shared, and combines these with the fecal shedding profile to estimate
# the effective number of shedders per day, to pair with number of ww detections.

# Note - this script will not run as we do not provide Hewitt et al's data, but
# we include it for complete transparency of our methods. The resultant fitted
# model parameters are provided for your use in downstream analyses or your own
# projects, and can be found in :
# "results/wastewater_sensitivity/fitted_params.csv"

# This script depends on:
### data from Hewitt et al (not provided in this repository)
### data/processed/shedding/fecal_shedding_profile_sc2.rds
### data/raw/sc2_shedding_profiles/shedding_He_plotdigitizer.csv
### R/data_analysis_functions/reshape_hewitt_data.R

# it produces:
### my_path/sc2_data_fecal_update.csv

# load libraries
source(here("R", "data_analysis_functions", "reshape_hewitt_data.R"))
library(tidyverse)
library(lubridate)
library(EpiSewer)
library(readxl)
library(rstan)

source(here("set_paths.R"))

################################################################################
#
################################################################################

## 1st: read in linelist of onset dates shared by authors
df1 <- readxl::read_xlsx(
  paste0(
    my_path, "/oct_2025_data_request.xlsx"
    )
  ) %>% mutate(
    day0 = if_else(asymptomatic == "No", # date of "day0" for each case
                        onset_date, lab_confirm_date)
    )%>% mutate(
      day0 = as.Date(day0)
      ) # POSIXIt --> date for easier manipulation later


## 2nd: read in data shared in SM Table S5
df2 <- readxl::read_xlsx(
  paste0(
    my_path, "/hewittSMtables.xlsx"),
  sheet = 5, range = "A1:I231"
  )%>% mutate(
    date = as.Date(date)
    )%>% # POSIXIt --> date for easier manipulation later
  pivot_wider(
    values_from = c(ww_detection, ww_log_genome_copies,
                  ww_av_ct, ww_positive_reps),
    names_from = ww_location
    )

## Load shedding distributions
shedding_dist_nasal <- read.csv(
  here(
    "data", "raw", "sc2_shedding_profiles", "shedding_He_plotdigitizer.csv")
)
names(shedding_dist_nasal) <- c("rel_day", "prob")

shedding_dist_fecal <- readRDS(
  here(
   "data","processed","shedding","fecal_shedding_profile_sc2.rds"
    )
)
names(shedding_dist_fecal) <- c("prob","rel_day")

## Reshape to get one row per day with number of shedders and effective shedders
daily_shedders_nasal <- compute_daily_shedding(df1,
                                         shedding_dist = shedding_dist_nasal,
                                         only_after_isolation = TRUE)
daily_shedders_fecal <- compute_daily_shedding(df1,
                                               shedding_dist = shedding_dist_fecal,
                                               only_after_isolation = TRUE)

## Attach daily shedders to wastewater detection data from SM
dat <- full_join(
  daily_shedders_nasal, df2,
  by = c("shedding_date" = "date")
  )
dat <- full_join(
  daily_shedders_fecal, dat,
  by = "shedding_date"
  )
dat <- dat %>%
  rename(
    "n_shedding.fecal" = `n_shedding.x`,
    "effective_shedding.fecal" = `effective_shedding.x`,
    "n_shedding.nasal" = `n_shedding.y`,
    "effective_shedding.nasal" = `effective_shedding.y`
    )

## do not keep dates where wastewater detection was not done
dat <- dat %>%
  filter(!is.na(ww_detection_WWTP | ww_detection_MIQF))

## Save
write.csv(dat, file = here(my_path,
                           "sc2_data_fecal_update.csv")
  )


