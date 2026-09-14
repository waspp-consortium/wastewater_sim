# This script simulates 1000 outbreaks for each of an array of parameter
# combinations defining viral characteristics, including viral load shed to
# wastewater, proportion of infections symptomatic, & probability of diagnosis.

# it depends on:
### results/estimates/ww_detection_params.rds - from script 4
### data/processes/shedding/fecal_shedding_profile_sc2.rds

# it produces:
### results/simulations/fig1/fig1_parx1-x2.rds - needed for Fig2A

# load libraries
library(wastewatchR)
library(tidyverse)
library(reshape2)
library(here)

################################################################################
#
################################################################################

# specify parameter grid
symptoms <- c(0.01, 0.05, seq(0.1, 1, by = 0.10))
shed_rel_sc2 <- 10^c(-6, -4, -2, -1, 0, 1, 2, 3)
delay_diagnosis_shape <- c(6)
prob_diagnosis <- c(0.1, 0.01, 1)

pars <- expand.grid(
  symptoms = symptoms,
  shed_rel_sc2 = shed_rel_sc2,
  delay_diagnosis_shape = delay_diagnosis_shape,
  prob_diagnosis = prob_diagnosis
  )
# add unique row id
pars <- pars %>% mutate(pars_row = row_number())

# read in detection curve parameters
# produced by scripts/4_refit_sensitivity_model_sc2_NZ.R
llpars <- as.data.frame(
  readRDS(
    here("results/estimates/ww_detection_params.rds")))$mean

# specify detection model
detection_params <- list(population = 10000, #from Hewitt et al.
                         duration = NA,
                         logistic_beta_0 = llpars[1],
                         logistic_beta_1 = llpars[2],
                         limit_of_detection = 0.01,
                         seed = NA)

# specify temporal shedding profile
# the shedding shape same as SARS-CoV-2 but scaled to reflect
# virus specific shedding amounts relative to sars-cov-2
shedding_dist <- readRDS(here(
  "data/processes/shedding/fecal_shedding_profile_sc2.rds"))$shedding

# generate n=R outbreaks for each parameter combination
R <- 1000  # number of stochastic replicates per parameter set

# seeds
set.seed(1234)  # master seed
base_seed <- 1234

# specify chunk of scenarios you want to run
## (best done in chunks due to run-time)
pars <- pars[249:288,]

results_df <- pars %>%
  pmap_dfr(function(symptoms, shed_rel_sc2, delay_diagnosis_shape, prob_diagnosis,pars_row) {

    res <- map_dfr(1:R, function(r) {

      # For each parameter row, run R stochastic outbreaks
      # Update seed so random but still reproducible using base seed
      seed_iter <- base_seed + pars_row * 1000 + r
      detection_params$seed <- seed_iter

      outbreak <- sim_single_outbreak(
        spillover_day = 1,
        index_case_ID = 1,
        mn_offspring = 0.99, # R0 just sub-critical
        disp_offspring = 1,
        max_gen = 100,
        index_cases = 1, # seed each outbreak with a single case
        generation_time_dist = function(n) rgamma(
          n, shape = 12, rate = 2
          ), # Table S6
        prob_symptomatic = symptoms,
        infection_to_onset_dist = function(n) rgamma(
          n, shape = 5.8, rate = 0.95
          ), # Table S6
        prob_severe = 1, # chose to aggregate prob severity, healthcare seeking and diagnosis --> prob of diagnosis
        prob_seek_healthcare_non_severe = 1, # chose to aggregate prob severity, healthcare seeking and diagnosis --> prob of diagnosis
        prob_seek_healthcare_severe = 1, # chose to aggregate prob severity, healthcare seeking and diagnosis --> prob of diagnosis
        onset_to_healthcare_dist = function(n) rgamma(
          n, shape = 14, rate = 2
          ),# Table S6
        prob_diagnosis = prob_diagnosis,
        healthcare_to_diagnosis_dist = function(n) rgamma(
          n, shape = delay_diagnosis_shape, rate = 2
          ), # Table S6
        initial_immune = 0
      )

      nts <- generate_number_shedding_time_series(
        outbreak,
        method = "effective_n_shedders",
        shedding_dist = shedding_dist,
        shedding_relative_SC2 = shed_rel_sc2
      )

      det <- calculate_wastewater_ttd(
        wastewater_number_shedding_time_series = nts,
        sampling_frequency = 7,
        sampling_method = "grab",
        detection_approach = "logistic_curve",
        detection_params
      )

      tibble(
        size = nrow(outbreak),
        day1WW = unlist(det$ttd),
        day1Clinic = if_else(
          sum(is.na(outbreak$time_diagnosis)) == nrow(outbreak),
          NA_real_,
          min(outbreak$time_diagnosis, na.rm = TRUE)
        ),
        total_symp = sum(outbreak$symptomatic),
        total_diag = sum(!is.na(outbreak$time_diagnosis)),
        iteration = r,
        pars_row = pars_row,
        symptoms = symptoms,
        prob_diagnosis = prob_diagnosis,
        delay_diagnosis_shape = delay_diagnosis_shape,
        shed_rel_sc2 = shed_rel_sc2,
        seed_used = seed_iter
      )

    })  # end map_dfr over R replicates

    # Print progress each time a full row of pars completes
    message("Completed pars_row: ", pars_row)

    res

  })

saveRDS(results_df, file = here("results/simulations/fig1/fig1_par249-288.rds"))
