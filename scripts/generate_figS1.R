# This script plots the numbers of infections per year for the different
# zoonotic archetypes

# It depends on:
# - results/simfig_s1.png generated whilst generating fig4C, S2 and S3

# It produces:
# - figures/figS1.png

# Load libraries
library(here)
library(tidyverse)

#-------------------------------------------------------------------------------

# Load model outputs
final_summary_metrics <- readRDS("results/simfig_s1.png")

## Plot
p_n_inf <- ggplot(data = final_summary_metrics%>%
                    filter(
                      freq == 7,
                      shedding_profile_name == "baseline_fecal")
                  )+
  geom_boxplot(aes(
    x = scenario, y = infections_day_365),
    fill = "grey80"
    )+
  theme_minimal()+
  xlab("")+
  ylab("Number of infections in year 1")+
  scale_x_discrete(labels = c(
    "Zoonotic archetype 1:\nhigh spillover rate\nlow R0",
    "Zoonotic archetype 2:\nepizootic with spillover\nno h-2-h transmission",
    "Zoonotic archetype 3:\nintense spike spillover risk\nh-2-h R0~1")
    )+
  coord_flip()+
  theme(text = element_text(size = 14))

ggsave(p_n_inf, filename = here("figures/figS1.png"),
       height = 4, width = 9, dpi = 300, unit = "in")
