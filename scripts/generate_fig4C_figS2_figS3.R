# This script generates Figure 4C, S2, and S3.

# It depends on:
# -  results/simulations/fig4C_simulations.rds

# It produces:
# - figures/fig4C.png
# - figures/figS2.png
# - figures/figS3.png

# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)
library(purrr)
library(here)

#-------------------------------------------------------------------------------

## load simulation results
final_summary_metrics <- readRDS(here(
  "results/simulations/fig4C_simulations.rds")
)

final_summary_metrics <- final_summary_metrics %>%
  mutate(scenario = factor(scenario, levels = c("lasv", "niv", "ebov")),
         freq = factor(freq, levels = c(1,7,14,28)),
         shedding_profile_name = names(shedding_profile),
         detected = case_when(
           !is.na(first_day_ww) & !is.na(first_day_cl) ~ "both",
           !is.na(first_day_ww) & is.na(first_day_cl) ~ "ww",
           is.na(first_day_ww) & !is.na(first_day_cl) ~ "cl",
           is.na(first_day_ww) & is.na(first_day_cl) ~ "neither",
           TRUE ~ as.character(NA)))%>%
  mutate(detected = factor(
    detected, levels = rev(c("ww", "both", "cl", "neither")))
    )

# split "both" into detected first in clinical and detected first in wastewater:
final_summary_metrics <- final_summary_metrics%>%
  mutate(detect_v2 = case_when(
    detected != "both" ~ detected,
    ((detected == "both") & (first_day_ww<=first_day_cl)) ~ "ww1st",
    ((detected == "both") & (first_day_cl<first_day_ww)) ~ "cl1st",
    TRUE ~ as.character(NA)),
    detect_v3 = case_when(detected != "both" ~ detected,
                          ((detected == "both") & (n_infections_ww<=n_infections_cl)) ~ "ww1st",
                          ((detected == "both") & (n_infections_cl<n_infections_ww)) ~ "cl1st",
                          TRUE ~ as.character(NA)))%>%
  mutate(detect_v2 = factor(detect_v2,
                            levels = rev(
                              c("ww", "ww1st", "cl1st", "cl", "neither")
                              )),
         detect_v3 = factor(detect_v3,
                            levels = rev(
                              c("ww", "ww1st", "cl1st", "cl", "neither")
                              )))
new_labels <- as_labeller(c(
  "lasv" = "Zoonotic archetype 1",
  "niv" = "Zoonotic archetype 2",
  "ebov" = "Zoonotic archetype 3",
  "0.01" = "1% probability of clinical diagnosis",
  "0.1" = "10% probability of clinical diagnosis",
  "0.001" = "0.1% probability of clinical diagnosis"
))

saveRDS(final_summary_metrics, file = here("results/simfig_s1.png"))

# Plot Figure 4C (performance according to infections prior to detection)
pS3 <- ggplot(final_summary_metrics%>%
                       filter(shedding_profile_name == "baseline_fecal"))+
  geom_bar(aes(x = freq, fill = detect_v2),
           position = "fill", color = "black")+
  scale_fill_manual(values = c( "salmon4",
                                "salmon3",
                                "#A6B8C8",
                                "#5B7C99",
                                "#D3D3D3" ),
                    labels = c("Wastewater only", "Wastewater first", "Clinical first", "Clinical only", "Neither"),
                    breaks = c("ww", "ww1st", "cl1st", "cl", "neither"),
                    name = "Detected by:")+
  facet_grid(prob_diagnosis~scenario, labeller = new_labels)+
  xlab("Frequency of wastewater sampling")+
  ylab("Proportion of model runs in which transmission was detected")+
  scale_x_discrete(labels = c("Daily", "Weekly", "Bi-weekly", "Monthly"))+
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 14),
        panel.grid = element_blank())
ggsave(pS3, filename = here("figures/figS3.png"),
       width = 9, height = 9, unit = "in")

## Plot 4C - just prob diag of 1%
p4C <- ggplot(final_summary_metrics%>%
                        filter(shedding_profile_name == "baseline_fecal",
                               prob_diagnosis == 0.01))+
  geom_bar(aes(x = freq, fill = detect_v2),
           position = "fill", color = "black", alpha=0.8)+
  scale_fill_manual(values = c( "salmon4",
                                "salmon3",
                                "#A6B8C8",
                                "#5B7C99",
                                "#D3D3D3" ),
                    labels = c("Wastewater only", "Wastewater first", "Clinical first", "Clinical only", "Neither"),
                    breaks = c("ww", "ww1st", "cl1st", "cl", "neither"),
                    name = "Detected by:")+
  facet_grid(prob_diagnosis~scenario, labeller = new_labels)+
  xlab("Frequency of wastewater sampling")+
  ylab("Proportion of model runs in which transmission was detected")+
  scale_x_discrete(labels = c("Daily", "Weekly", "Bi-weekly", "Monthly"))+
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 14),
        panel.grid = element_blank())
ggsave(p4C, filename = here("figures/fig4C.png"),
       width = 9, height = 4, unit = "in")

# Plot S2 - effect of shedding profile
pS2 <- ggplot(final_summary_metrics)+
  geom_bar(aes(x = shedding_profile_name, fill = detect_v2),
           position = "fill", color = "black")+
  scale_fill_manual(values = c( "salmon4",
                                "salmon3",
                                "#A6B8C8",
                                "#5B7C99",
                                "#D3D3D3" ),
                    labels = c("Wastewater only", "Wastewater first", "Clinical first", "Clinical only", "Neither"),
                    breaks = c("ww", "ww1st", "cl1st", "cl", "neither"),
                    name = "Detected by:")+
  facet_grid(prob_diagnosis~scenario, labeller = new_labels)+
  xlab("Shedding profile shape")+
  scale_x_discrete(labels = c("baseline\nfecal",
                              "spikier\nfecal",
                              "baseline\nnasal",
                              "spikier\nnasal",
                              "flat"))+
  ylab("Proportion of model runs in which transmission was detected")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 14))
ggsave(pS2, filename = here("figures/figS2.png"),
       width = 10, height = 10, unit = "in")

