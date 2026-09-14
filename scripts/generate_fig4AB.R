# This script plots figures 4A and 4B

# It depends on the outputs of "simulations_fig4AB.R":
# - results/simulations/fig4B_summary_prob_diag001.rds
# - results/simulations/fig4B_summary_prob_diag0001.rds
# - results/simulations/fig4B_summary_prob_diag01.rds
# - outbreaks for 4A incidence: "results/simulations/fig4B_probdiag01.rds"

# It produced:
# - figures/fig4A.png
# - figures/fig4B.png

# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)
library(purrr)
library(here)

# ------------------------------------------------------------------------------

## bring together different runs
final_summary_metrics001 <- readRDS(here(
  "results/simulations/fig4B_summary_prob_diag001.rds")
)
final_summary_metrics0001 <- readRDS(here(
  "results/simulations/fig4B_summary_prob_diag0001.rds")
)
final_summary_metrics01 <- readRDS(here(
  "results/simulations/fig4B_summary_prob_diag01.rds")
)

dat <- bind_rows(final_summary_metrics0001,
                 final_summary_metrics001,
                 final_summary_metrics01)

pdat <- dat %>%
  dplyr::mutate(detected = if_else(!is.na(first_day_ww) | !is.na(first_day_cl), "Y", "N"),
                WW_only_or_early_days = case_when(
                  !is.na(first_day_ww) & is.na(first_day_cl) ~ "Y",
                  (!is.na(first_day_ww) & !is.na(first_day_cl) & first_day_ww<first_day_cl) ~ "Y",
                  is.na(first_day_ww) & !is.na(first_day_cl) ~ "N",
                  (!is.na(first_day_ww) & !is.na(first_day_cl) & first_day_ww>=first_day_cl) ~ "N",
                  is.na(first_day_ww) & is.na(first_day_cl) ~ as.character(NA),
                  TRUE ~ as.character(NA)),
                WW_only_or_early_infs = case_when(
                  !is.na(n_infections_ww) & is.na(n_infections_cl) ~ "Y",
                  (!is.na(n_infections_ww) & !is.na(n_infections_cl) & n_infections_ww<n_infections_cl) ~ "Y",
                  is.na(n_infections_ww) & !is.na(n_infections_cl) ~ "N",
                  (!is.na(n_infections_ww) & !is.na(n_infections_cl) & n_infections_ww>=n_infections_cl) ~ "N",
                  is.na(n_infections_ww) & is.na(n_infections_cl) ~ as.character(NA),
                  TRUE ~ as.character(NA)))%>%
  group_by(
    scenario, prob_symptomatic, shedding_relative_SC2, prob_diagnosis, freq
    )%>%
  dplyr::summarise(n = n(),
                   n_detected = sum(detected == "Y"),
                   percentage_detected = sum(detected == "Y")/n,
                   percentage_detected_WW = 100 * sum(!is.na(first_day_ww))/n,
                   percentage_detected_clinical = 100 * sum(!is.na(first_day_cl))/n,
                   percentage_of_detected_WW_best_days = 100 * sum(
                     WW_only_or_early_days == "Y", na.rm = T
                     )/sum(
                       detected == "Y", na.rm= T
                       ),
                   percentage_of_detected_WW_best_infs = 100 * sum(
                     WW_only_or_early_infs == "Y", na.rm = T
                     )/sum(
                       detected == "Y", na.rm= T)
                   )%>%
  dplyr::mutate(WW_best_days = if_else(
    percentage_of_detected_WW_best_days>50, "Y", "N"),
                WW_best_infs = if_else
    (percentage_of_detected_WW_best_infs>50, "Y", "N")
    )


# extract points at which WW_best_days flips to "Y"
thresholds <- pdat %>%
  dplyr::group_by(shedding_relative_SC2,
                  prob_diagnosis, scenario) %>%
  dplyr::summarise(
    threshold_symptoms = max(
      ifelse(WW_best_days == "Y", prob_symptomatic, 0),
      na.rm = TRUE
    ),
    .groups = "drop"
  )%>%
  mutate(scenario = factor(scenario, levels = c("lasv",
                                                "niv",
                                                "ebov")))

#-------------------------------------------------------------------------------
# plot Figure 4B

p4B <- ggplot(
  thresholds,
  aes(
    x = shedding_relative_SC2,
    y = 100*threshold_symptoms,
    fill = scenario,
    colour = scenario,
    group = prob_diagnosis
  )
) +
  geom_ribbon(
    aes(ymin = 0, ymax = 100*threshold_symptoms),
    alpha = 0.5, col = NULL
  ) +
  #geom_line() +
  scale_fill_manual(values = c(
    "#5C8A5C",
    "#8A6FA8",
    "#B89B4A"
  ))+
  facet_wrap(~scenario, nrow = 1) +
  theme_classic()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  ylab("Infections with clinical symptoms (%)")+
  xlab("Mean viral load shed (copies/day) relative to SARS-CoV-2")+
  theme(text = element_text(size = 14),
        legend.position = "none")+
  scale_x_log10(
    limits = c(1e-6, 1e+3),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
ggsave(p4B, filename = here(
  "figures/fig4B.png"),
  dpi = 300, width = 9, height = 4, unit = "in")

#-------------------------------------------------------------------------------
# Plot figure 4A

# load weekly incidence (can pick any batch as prob diagnosis does not effect
# incidence)
outbreaks <- readRDS(file = here(
  "results/simulations/fig4B_probdiag01.rds")
)
replicate_tbl <- outbreaks %>%
  unnest(sim) %>%
  group_by(scenario, prob_diagnosis, prob_symptomatic, replicate) %>%
  nest()

max_time <- 3 * 365  # 3 years in days

# extract weekly incidence
weekly_incidence <- replicate_tbl %>%

  # 0. filter so prob_symptomatic == 1 (just to cut down time)
  filter(prob_symptomatic == 1) %>%

  # 1. unpack nested simulation data
  unnest(data) %>%

  # 2. create weekly bins
  mutate(week = floor(time_infection / 7) * 7) %>%

  # 3. count infections per week per simulation setting
  count(
    scenario,
    replicate,
    prob_diagnosis,
    prob_symptomatic,
    week,
    name = "incidence"
  ) %>%

  # 4. enforce full weekly time grid (0–3 years)
  group_by(scenario, replicate, prob_diagnosis, prob_symptomatic) %>%
  complete(
    week = seq(0, max_time, by = 7),
    fill = list(incidence = 0)
  ) %>%
  ungroup()

scale_factor <- 333 #(double axes)

weekly_incidence <-
  weekly_incidence %>%
  mutate(scenario = factor(scenario,
                           levels = c(
                             "lasv", "niv", "ebov"))
         )
p4A <- ggplot() +
  # incidence (left axis)
  geom_line(
    data = weekly_incidence,
    aes(
      x = week,
      y = incidence,
      group = interaction(replicate, prob_diagnosis, prob_symptomatic)
    ),
    alpha = 0.1
  ) +

  # spillover (right axis, scaled)
  geom_line(
    data = spill,
    aes(
      x = t,
      y = spillover_rate * scale_factor,
      colour = scenario
    ),
    linewidth = 1
  ) +
  facet_wrap(~ scenario, ncol = 3) +
  scale_colour_manual(values = c(
    lasv = "#5C8A5C",
    niv  = "#8A6FA8",
    ebov = "#B89B4A"
  )) +
  scale_y_continuous(
    name = "Incidence",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Spillover rate"
    )
  ) +
  labs(x = "Day") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 365))+
  theme(legend.position = "none",
        text = element_text(size = 14),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(p4A, filename = here(
  "figures/fig4A.png"
  ), dpi = 300, width = 9, height = 2, unit = "in")



