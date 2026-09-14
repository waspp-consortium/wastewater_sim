# This script prepares alternative shedding profiles with same total viral load
# but different distribution over time. It's used to assess robustness of
# results to assumed shedding profiles (shown in Figure S2).

# It depends on:
# -  SARS-CoV-2 profile: data/processed/shedding/fecal_shedding_profile_params_sc2.rds

# It produces:
# - List of profiles for sensitivity analysis: data/processed/shedding/all_profiles.rds
# - Normalised AUC for SARS-CoV-2 fecal shedding profile: results/AUC_wannigama.rds

# Load libraries
library(EpiSewer)
library(tidyverse)

#-------------------------------------------------------------------------------

# Load baseline params (sC2 fecal shedding profile)
wannigama_params <- readRDS(here(
  "data/processed/shedding/fecal_shedding_profile_params_sc2.rds")
  )
# Compute mean and sd for ease of reshaping later
gamma_mean <- wannigama_params[["k"]] * wannigama_params[["theta"]]
gamma_sd   <- wannigama_params[["theta"]] * sqrt(wannigama_params[["k"]])

# Define the day range
days <- -3:21
shift <- -min(days)
max_day <- max(days + shift)

# Helper function - do not normalise yet
generate_raw <- function(gamma_mean, gamma_sd, days, max_day) {
  shedding <- EpiSewer::get_discrete_gamma(
    gamma_mean = gamma_mean,
    gamma_sd   = gamma_sd,
    max        = max_day
  )

  data.frame(day = days, shedding = shedding)
}

# Create raw distributions
baseline_raw <- generate_raw(gamma_mean, gamma_sd, days, max_day)
spikier_raw <- generate_raw(gamma_mean, gamma_sd * 0.25, days, max_day)

# Normalise to ensure peak baseline = 1
baseline <- baseline_raw %>%
  mutate(shedding = shedding / max(shedding))
target_auc <- sum(baseline$shedding)

spikier <- spikier_raw %>%
  mutate(shedding = shedding / sum(shedding) * target_auc)

# Prepare the same using a shedding profile with a much earlier peak
# (similar to nasal shedding in SC2)
nasal_params <- readRDS(here(
  "data/processed/shedding/nasal_shedding_profile_params_sc2.rds")
)

# Compute mean and sd for ease of reshaping later
gamma_mean <- nasal_params[["k"]] * nasal_params[["theta"]]
gamma_sd   <- nasal_params[["theta"]] * sqrt(nasal_params[["k"]])

# Define the day range
days <- -3:21
shift <- -min(days)
max_day <- max(days + shift)

# Create raw distributions
baseline_raw <- generate_raw(gamma_mean, gamma_sd, days, max_day)
spikier_raw <- generate_raw(gamma_mean, gamma_sd * 0.25, days, max_day)

# Normalise to ensure peak baseline = 1
baseline_nasal <- baseline_raw %>%
  mutate(shedding = shedding / max(shedding))
target_auc_nasal <- sum(baseline$shedding)

spikier_nasal <- spikier_raw %>%
  mutate(shedding = shedding / sum(shedding) * target_auc_nasal)

# also prepare flat extreme
flat_raw <- baseline_raw %>% mutate(shedding = 1)
flat <- flat_raw %>%
  mutate(shedding = shedding / sum(shedding) * target_auc_nasal)

# bring together
all_profiles <- bind_rows(
  baseline %>% mutate(profile = "Baseline"),
  spikier  %>% mutate(profile = "Spikier"),
  baseline_nasal %>% mutate(profile = "Baseline nasal"),
  flat  %>% mutate(profile = "Flat"),
  spikier_nasal  %>% mutate(profile = "Spikier nasal")
)

# plot as a sanity check
ggplot(all_profiles, aes(x = day, y = shedding, colour = profile)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Comparison of Viral Shedding Profiles",
    subtitle = "Baseline peaks at 1; all profiles have equal total shedding",
    x = "Days PSO",
    y = "Normalised viral shedding",
    colour = "Profile"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  scale_x_continuous(breaks = seq(min(all_profiles$day),
                                  max(all_profiles$day), by = 2))

# Save all alternative profiles
list_profiles <- list(baseline_fecal = baseline,
                      spikier_fecal = spikier,
                      baseline_early = baseline_nasal,
                      flat = flat,
                      spikier_early = spikier_nasal)
saveRDS(list_profiles, here(
  "data/processed/shedding/all_profiles.rds")
)

# Save AUC for anchoring shedding estimates
saveRDS(file = here("results/AUC_wannigama.rds"), target_auc)
