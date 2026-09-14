# This script runs the simulations needed to plot Figure 4C S2 amd S3 - looking
# at how viral transmission dynamics, summarised as 3 different zoonotic
# archetypes effect the value added by wastewater surveillance.

# Unlike the simulations for 4AB, multiple frequencies of wastewater sampling
# and multiple shedding profile shapes are used, but only 1 baseline set of
# viral characteristics (shedding amounts and clinical symptom proportion)
# based on a SARS-CoV-2-like pathogen.

# It depends on:
# - R/simulation_functions/estimate_clinical_detection_prob.R

# It produces:
# - results/simulations/fig4C_simulations.rds

# load libraries
library(wastewatchR)
library(tidyverse)
library(reshape2)
library(igraph)
library(ggraph)
library(cowplot)
library(purrr)
source("R/simulation_functions/estimate_clinical_detection_prob.R")

#-------------------------------------------------------------------------------

# define the zoonotic archetypes
lasv_pars <- c("mean_sr" = 0.3, "P" = 365,
               "tmax" = 183, "d" = 3, "sigma" = 63,
               "b" = solve_b(x=0.3,tmax=183,sigma=60,d=3),
               "R0" = 0.2, "disp_offspring" = 2)
niv_pars <- c("P" = 5*365,
              "tmax" = 183, "d" = 300, "sigma" = 35,
              "b" = 0.001,
              "R0" = 0, "disp_offspring" = 1)
ebov_pars <- c("P" = 5*365,
               "tmax" = 183, "d" = 240, "sigma" = 3.5,
               "b" = 0.001,
               "R0" = 0.95, "disp_offspring" = 2)
scenarios <- list(
  lasv = lasv_pars,
  niv  = niv_pars,
  ebov = ebov_pars
)

# varying probability of diagnosis
prob_diagnosis <- c(0.01, 0.1, 0.001)

param_grid <- tidyr::crossing(
  scenario = names(scenarios),
  prob_diagnosis = prob_diagnosis
)

# set global params for clinical delays and shedding
standard_pars <- c(prob_symptomatic = 1,
                   prob_seek_healthcare = 1,
                   generation_time_dist = function(n)
                   { rgamma(n, shape = 12, rate = 2) },
                   infection_to_onset_dist = function(n)
                   { rgamma(n, shape = 6, rate = 2) },
                   prob_severe = 1,
                   prob_seek_healthcare_non_severe = 0,
                   onset_to_healthcare_dist = function(n)
                   { rgamma(n, shape = 9, rate = 1.5) },
                   healthcare_to_diagnosis_dist = function(n)
                   { rgamma(n, shape = 6, rate = 2) },
                   initial_immune = 0 # CHECK THIS
)

# run spillover models
results <- imap_dfr(
  scenarios,
  function(pars, scenario_name){

    map_dfr(
      1:1000,
      function(rep){

        spillover(
          time = 3*365,
          specify = "spillover_rate",
          seasonal_period = pars["P"],
          spillover_rate_pars = list(
            tmax = pars["tmax"],
            b = pars["b"],
            d = pars["d"],
            sigma = pars["sigma"]
          ),
          prevalence_pars = NULL,
          contact_rate_pars = NULL,
          p_infection = NULL
        ) %>%
          mutate(
            scenario = scenario_name,
            replicate = rep
          )

      }
    )

  }
)

# run h-2-h transmission models

### one row per spillover event
spillover_events <- results %>%
  filter(spillovers > 0) %>%
  uncount(spillovers) %>%
  group_by(scenario, replicate) %>%
  mutate(event_id = row_number()) %>%
  ungroup()
spillover_events <- spillover_events %>% # join w/ param_grid for sens. analysis
  left_join(param_grid, by = "scenario")

outbreaks <- spillover_events %>%
  mutate(
    sim = pmap(
      list(t, event_id, scenario, prob_diagnosis),
      function(spillover_day, index_case_ID, sc, prob_diag){

        pars <- scenarios[[sc]]

        sim_single_outbreak(
          spillover_day = spillover_day,
          index_case_ID = index_case_ID,
          mn_offspring = pars["R0"],
          disp_offspring = pars["disp_offspring"],
          max_gen = Inf,
          index_cases = 1,
          generation_time_dist = standard_pars[["generation_time_dist"]],
          prob_symptomatic = standard_pars[["prob_symptomatic"]],
          infection_to_onset_dist = standard_pars[["infection_to_onset_dist"]],
          prob_severe = 1,
          prob_seek_healthcare_non_severe = 0,
          prob_seek_healthcare_severe = standard_pars[["prob_seek_healthcare"]],
          onset_to_healthcare_dist = standard_pars[["onset_to_healthcare_dist"]],
          prob_diagnosis = prob_diag,
          healthcare_to_diagnosis_dist = standard_pars[["healthcare_to_diagnosis_dist"]],
          initial_immune = 0,
          seed = sample.int(1e6, 1)
        )
      },
      .progress =TRUE
    )
  )

# Collapse into linelist per replicate per scenario
replicate_tbl <- outbreaks %>%
  unnest(sim) %>%
  group_by(scenario, prob_diagnosis, replicate) %>%
  nest()

# Set surveillance parameters:
logistic_params <- readRDS(here(
  "results/estimates/ww_detection_params.rds")
)
frequency_sampled <- c(1, 7, 14, 28)

# Shedding profiles
shedding_dist_list <- readRDS(
  "data/processed/shedding/all_profiles.rds")

shedding_pars <- tidyr::crossing(
  shedding_profile = shedding_dist_list
)
sampling_pars <- tidyr::crossing(
  freq = frequency_sampled
)
replicate_tbl <- replicate_tbl %>%
  tidyr::crossing(shedding_pars)

# wastewater detection parameters
wastewater_params <- list(
  sampling_method = "grab",
  detection_approach = "logistic_curve",
  detection_params = list(
    population = 100000,
    logistic_beta_0 = logistic_params[1,1],
    logistic_beta_1 = logistic_params[2,1],
    limit_of_detection = 0.01,
    seed = 123
  )
)

# 1. Generate shedding time series
shedding_time_series <- replicate_tbl %>%
  mutate(
    shedding_ts = map2(
      data,
      shedding_profile,
      ~ generate_number_shedding_time_series(
        branching_process_output = .x,
        shedding_dist = .y$shedding,
        shedding_relative_SC2 = 1,
        method = "effective_n_shedders",
        duration = NA
      ) %>%
        ungroup() %>%
        arrange(day) %>%
        mutate(cumulative_infections = cumsum(new_infections))
    )
  )

shedding_time_series <- shedding_time_series %>%
  tidyr::crossing(sampling_pars) # expand over frequency options

# 2. Wastewater detection
wastewater_detection_results <- shedding_time_series %>%
  ungroup() %>%
  mutate(
    ww_detection = map2(
      shedding_ts,
      freq,
      ~ calculate_wastewater_ttd(
        wastewater_number_shedding_time_series = .x,
        sampling_frequency = .y,
        sampling_method = wastewater_params$sampling_method,
        detection_approach = wastewater_params$detection_approach,
        detection_params = wastewater_params$detection_params
      )
    ),
    ttd = map(ww_detection, "ttd"))

# 3. Summarise detection metrics
final_summary_metrics <- wastewater_detection_results %>%
  mutate(

    # 1. first day detected in ww from ttd
    first_day_ww = map_dbl(ttd, ~ .x[[1]]),

    # 2. first day detected clinically = min time_diagnosis
    first_day_cl = map_dbl(data, ~ {
      x <- .x$time_diagnosis
      if (all(is.na(x))) return(NA_real_)
      min(x, na.rm = TRUE)
    }),

    # 3. First day of human infection (min(t) where et indicates infection)
    first_day_infection = map_dbl(data, ~ {
      infection_times <- .x %>%
        pull(time_infection)

      if (length(infection_times) == 0 || all(is.na(infection_times))) {
        NA_real_
      } else {
        min(infection_times, na.rm = TRUE)
      }
    }),

    # 4. infections at detection time in ww
    n_infections_ww = pmap_dbl(
      list(shedding_ts, first_day_ww),
      function(df, t) {

        if (is.na(t)) return(NA_real_)

        val <- df %>%
          filter(day == t) %>%
          pull(cumulative_infections)

        if (length(val) == 0) return(NA_real_)

        first(val)
      }
    ),

    # 5. cl infections at first clinical detection day
    n_infections_cl = pmap_dbl(
      list(shedding_ts, first_day_cl),
      function(df, t) {

        if (is.na(t)) return(NA_real_)

        val <- df %>%
          filter(day == floor(t)) %>%
          pull(cumulative_infections)

        if (length(val) == 0) return(NA_real_)

        first(val)
      }
    ),
    # 6. Extract umulative infections at end of each year
    infections_day_365 = map_dbl(shedding_ts, ~ {
      idx <- which.min(abs(.x$day - 365))
      .x$cumulative_infections[idx]
    }),

    infections_day_730 = map_dbl(shedding_ts, ~ {
      idx <- which.min(abs(.x$day - 730))
      .x$cumulative_infections[idx]
    }),

    infections_day_1095 = map_dbl(shedding_ts, ~ {
      idx <- which.min(abs(.x$day - 1095))
      .x$cumulative_infections[idx]
    })
  ) %>%
  select(scenario, prob_diagnosis, freq, shedding_profile,
         replicate, first_day_infection,
         first_day_ww, first_day_cl,
         n_infections_ww, n_infections_cl,
         infections_day_365, infections_day_730, infections_day_1095)

# save
saveRDS(final_summary_metrics, file = "results/simulations/fig4C_simulations.rds")
