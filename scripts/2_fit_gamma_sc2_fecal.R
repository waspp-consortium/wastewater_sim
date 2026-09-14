# This script processes the fecal shedding profile data extracted from
# Wannigama et al and fits a gamma distribution for use in downstream modelling

# it depends on:
### - data/raw/sc2_shedding_profiles/wannigama_fecal_shedding_PSO.csv

# it produces:
### data/processed/shedding/fecal_shedding_profile_sc2.rds
### data/processed/shedding/fecal_shedding_profile_params_sc2.rds"

# load libraries
library(EpiSewer) # has a nice function for discretising the gamma
library(tidyverse)

################################################################################
#
################################################################################

# Gamma function, where A is amplitude
gamma_fun <- function(day, k, theta, A) {
  A * dgamma(day, shape = k, scale = theta)
}

# Read in shedding profile based on fecal shedding from Wannigama et al.,
# https://doi.org/10.1016/S1473-3099(24)00155-5 extracted using plotdigitizer
shedding_wannigama <- read.csv(
  here(
    "data/raw/sc2_shedding_profiles/wannigama_fecal_shedding_PSO.csv"
    )
  )

# Get the mean viral load shed at each time point across all variants excluding
# JN.1 and BA.286 (higher shedding, not circulating during period relevant to us)
shedding_wannigama <- shedding_wannigama %>%
  filter(
    variant != "JN1" & variant != "BA286"
    )%>%
  group_by(x)%>%
  summarise(
    y = mean(y)
    )

# Shift to avoid x<=0 for gamma distribution
epsilon <- 1e-6
shift <- min(
  shedding_wannigama$x
  )
shedding_wannigama$day_shifted <-
  shedding_wannigama$x - shift + epsilon

# Fit gamma
fit <- nls(
  y ~ gamma_fun(day_shifted, k, theta, A),
  data = shedding_wannigama,
  start = list(
    k = 3,    # shape
    theta = 3,    # scale
    A = max(shedding_wannigama$y) * 10
  ),
  algorithm = "port", # more stable
  lower = list(k=0.1, theta=0.1, A=0)   # avoid invalid parameter combos
)

# Extract fitted parameters
params <- coef(fit)
shape <- params["k"]
scale <- params["theta"]

# Discretise and normalise
wannigama_dist <- data.frame(
  shedding = EpiSewer::get_discrete_gamma(
    gamma_mean = shape * scale,
    gamma_sd   = scale * sqrt(shape),
    max = max(shedding_wannigama$x) - shift),
  day = -3:21)
wannigama_dist <- wannigama_dist %>%
  mutate(
    shedding = shedding/max(wannigama_dist$shedding)
    )

# Save shedding profile parameters for use in Hewitt refits and modelling
saveRDS(file = here(
  "data/processed/shedding/fecal_shedding_profile_sc2.rds"), wannigama_dist
)
saveRDS(file = here(
  "data/processed/shedding/fecal_shedding_profile_params_sc2.rds"), params
  )

# Plot data and discretised fitted profile as a sense check
p <- ggplot(
  shedding_wannigama
  )+
  geom_bar(
    data = wannigama_dist,
    aes(
      x = day, y = shedding
      ),
    stat = "identity",
    fill = "grey80",
    )+
  geom_point(
    aes(
      x = x, y = y/(max(shedding_wannigama))
    )
  )+
  theme_minimal()+
  xlab("days PSO")+
  ylab("Normalised fecal shedding")

print(p)
