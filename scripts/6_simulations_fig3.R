# This script provides an example (used to generate Figure 3 in the manuscript)
# of how to use wastewatchR to simulate spillover, onward transmission, and
# surveillance for an emerging virus. In this case the virus has a high spillover
# rate and low R0 similar to LASV, with a shedding profile and clinical presentation
# and delays similar to SARS-CoV-2

# It depends on:
### Detection model parameters: results/estimates/ww_detection_params.rds
### Shedding profile: data/processes/shedding/fecal_shedding_profile_sc2.rds
### function "calculate_clinical_detection_prob" in R/simulation_functions

# It produces:
### Simulation to be plotted in generate_fig3:
##### results/simulations/fig3.rds
##### results/simulations/fig3_linelist.rds

# load libraries
library(wastewatchR)
library(tidyverse)
library(reshape2)
library(igraph)
library(ggraph)
library(cowplot)
library(here)

source(here("R/simulation_functions/estimate_clinical_detection_prob.R"))

# 1. Simulate spillovers over a 1-yr period ------------------------------------
# ------------------------------------------------------------------------------
# Here assume a seasonal spillover rate that varies with an annual periodicity
# peaking on day 183 of the year, with an average daily spillover rate of 0.3
# in our chosen catchment area of 100,000 people

catchment_size <- 100000

sp <- spillover(time = 365, specify = "spillover_rate",
                seasonal_period = 365,
                spillover_rate_pars = list("tmax" = 183,
                                           "b" = solve_b(
                                             x=0.3,tmax=183,sigma=60,d=6
                                             ),
                                           "d" = 3,
                                           "sigma" = 63))%>%
  filter(spillovers >0)

# account for possibility of multiple spillovers on the same day
sp <- sp %>%
  uncount(spillovers) %>%
  rowid_to_column("ID")

# 2. For each spillover, simulate onward transmission based aon assumed R0 -----
# ------------------------------------------------------------------------------

# List of spillover events that could initiate onward transmission
events <- vector("list", length = dim(sp)[1])

# Specify clinical surveillance parameters
## Collapse all probabilities into aggregate "probability of diagnosis"
prob_symptomatic <- 1
prob_seek_healthcare <- 1
prob_diagnosis <- 0.01

## delay infection-->onset
shape_onset <- 6
rate_onset <- 2

## delay onset--> healthcare
shape_hc <- 9
rate_hc <- 1.5

## delay healthcare-->diagnosis
shape_diag <- 6
rate_diag <- 2

# Simulate outbreaks following each spillover event
## Specify offspring distribution and delays:

for(s in 1:(length(events))){
  events[[s]] <- sim_single_outbreak(spillover_day = sp$t[s],
                                     index_case_ID = sp$ID[s],
                                     mn_offspring = 0.2, # low R0 of 0.2
                                     disp_offspring = 2, # >1 --> overdispersion
                                     max_gen = Inf,
                                     index_cases = 1,
                                     generation_time_dist = function(n)
                                     { rgamma(n, shape = 12, rate = 2) },
                                     prob_symptomatic = prob_symptomatic,
                                     infection_to_onset_dist = function(n)
                                     { rgamma(n, shape = shape_onset, rate = rate_onset) },
                                     prob_severe = 1, # given symptomatic
                                     prob_seek_healthcare_non_severe = 0,
                                     prob_seek_healthcare_severe = prob_seek_healthcare,
                                     onset_to_healthcare_dist = function(n)
                                     { rgamma(n, shape = shape_hc, rate = rate_hc) },
                                     prob_diagnosis = prob_diagnosis,
                                     healthcare_to_diagnosis_dist = function(n)
                                     { rgamma(n, shape = shape_diag, rate = rate_diag) },
                                     initial_immune = 0, # no prior immunity
                                     seed = 102)
}

# merge linelists
linelist <- bind_rows(events)%>%
  mutate(type = if_else(infector == "animal", "Animal-to-human", "Human-to-human"))


# 3. Wastewater detection ------------------------------------------------------
# ------------------------------------------------------------------------------

# Load your assumed shedding distribution here
## as a placeholder we used one based on fecal shedding profile of SARS-CoV-2

shedding_dist <- readRDS(here(
  "data/processes/shedding/fecal_shedding_profile_sc2.rds")
  )

# Generate time series of number of individuals shedding each day
## Choose whether the model is based on the total number of people shedding
## (i.e. anyone shedding counts as one shedder), or on the effective number
## of shedders, which accounts for changes in shedding amount over time since
## infection. Under the latter approach, a person on their peak day of shedding
## counts as 1 shedder, while someone shedding at 10% of their peak amount
## counts as 0.1 shedders.
nts <- generate_number_shedding_time_series(linelist,
                                            method = "effective_n_shedders",
                                            shedding_dist = shedding_dist$shedding,
                                            shedding_relative_SC2 = 1)

# Specify your chosen model describing the relationship between your measure of
# shedding (number of shedders or effective shedders) and probabilty of detection
## Here we use a model based on data from SARS-CoV-2 at a quarantine facility in
## New Zealand.
logistic_params <- readRDS(here(
  "results/estimates/ww_detection_params.rds")
)
detection_params <- list(
  population = catchment_size,
  logistic_beta_0 = logistic_params[1,1],
  logistic_beta_1 = logistic_params[2,1],
  limit_of_detection = 0.01, # Specify LOD in terms of effective shedders
  seed = 123
)

# simulate wastewater detection - eventually want to take this out the package and combine this and clinical prob functions?
det <- calculate_wastewater_ttd(wastewater_number_shedding_time_series = nts,
                                sampling_frequency = 7,
                                sampling_method = "grab",
                                detection_approach = "logistic_curve",
                                detection_params)

# compute probability of clinical detection each day based on delay distributions and asusmed probs

det_clinical <- calculate_clinical_detection_prob(nts = nts,
                                                  det = det,
                                                  delay_onset = dgamma(0:21, shape = 12, rate = 2),
                                                  delay_seek = dgamma(0:21, shape = 9, rate = 1.5),
                                                  delay_diag = dgamma(0:21, shape = 6, rate = 2),
                                                  prob_symptomatic = prob_symptomatic,
                                                  prob_seek_healthcare = prob_seek_healthcare,
                                                  prob_diagnosis = prob_diagnosis)

saveRDS(det_clinical, here("results/simulations/fig3.rds"))
saveRDS(linelist, here("results/simulations/fig3_linelist.rds"))
