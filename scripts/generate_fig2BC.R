# This script generates figure 2 panel B and C

# It depends on:
# - results of "5_simulations_fig2.R"

# It produces:
# - "figures/fig2BC.png"


# Load libraries
library(tidyverse)
library(cowplot)
library(purrr)
library(here)

#-------------------------------------------------------------------------------

## read in results (batches for convenience of run times)
r1 <- readRDS(
  here("results/simulations/fig1/fig1_par1-24.rds")
)
r2 <- readRDS(
  here("results/simulations/fig1/fig1_par25-36.rds")
)
r3 <- readRDS(
  here("results/simulations/fig1/fig1_par37-66.rds")
)
r4 <- readRDS(
  here("results/simulations/fig1/fig1_par67-96.rds")
)
r5 <- readRDS(
  here("results/simulations/fig1/fig1_par97-136.rds")
)
r6 <- readRDS(
  here("results/simulations/fig1/fig1_par137-166.rds")
)
r7 <- readRDS(
  here("results/simulations/fig1/fig1_par167-200.rds")
)
r8 <- readRDS(
  here("results/simulations/fig1/fig1_par201-230.rds")
)
r9 <- readRDS(
  here("results/simulations/fig1/fig1_par231-248.rds")
)
r10 <- readRDS(
  here("results/simulations/fig1/fig1_par249-288.rds")
)

## merge
rtot <- dplyr::bind_rows(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)

rtot <- rtot %>%
  mutate(symptoms = round(
    symptoms,2)) # fix floating point issue

## Compute detection types
rtot <- rtot %>%
  filter(size>=10, delay_diagnosis_shape == 6)%>%
  mutate(detection = case_when(!is.na(day1WW) & !is.na(day1Clinic) ~ "both",
                               !is.na(day1WW) & is.na(day1Clinic) ~ "ww",
                               is.na(day1WW) & !is.na(day1Clinic) ~ "cl",
                               is.na(day1WW) & is.na(day1Clinic) ~ "neither"),
         time_saved = day1Clinic - day1WW)

## Select three example scenarios
rtot <- rtot %>%
  filter((prob_diagnosis == 0.1 & (((shed_rel_sc2 == 1e+00) & (symptoms == 0.7)) |
                                   ((shed_rel_sc2 == 1e+03) & (symptoms == 0.01)) |
                                   ((shed_rel_sc2 == 1e-02) & (symptoms == 0.4)))))

## Make detection a factor ordered for plot clarity
rtot <- rtot%>%
  mutate(detection = factor(detection,
                            levels = c("neither", "cl", "both", "ww")),
         profile = paste0("shedding=", shed_rel_sc2, ", symp=", symptoms))

## Plot B
pb <- ggplot(rtot%>%filter(profile!="shedding=0.01, symp=0.2"), aes(x = factor(profile), fill = detection)) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "",
    y = "Outbreaks detected (%)",
    fill = "Method\ndetected"
  ) +
  scale_x_discrete(labels = c(#"Low shedding\nLow symptomatic\ne.g. DENV/CCHFV/LASV",
                              "Low shedding\nMid symptomatic\ne.g. ZIKV",
                              "Mid shedding\nHigh symptomatic\ne.g. SARS-CoV-2",
                              "High shedding\nLow symptomatic\ne.g. Poliovirus")) +
  scale_fill_manual(values = c( "salmon4",
                                "#5B7C99",
                                "grey30",
                                "#D3D3D3" ),
                    labels = c("Wastewater", "Clinical", "Both", "Neither"),
                    breaks = c("ww", "cl", "both", "neither"))+
  theme_classic()+
  theme(legend.position = "top",
        text = element_text(size = 14))


## Plot C
custom_labels <- c("shedding=0.01, symp=0.2" = "Low shedding\nLow symptomatic\ne.g. DENV/CCHFV/LASV",
                   "shedding=0.01, symp=0.4"  = "Low shedding\nMid symptomatic\ne.g. ZIKV",
                   "shedding=1, symp=0.7" = "Mid shedding\nHigh symptomatic\ne.g. SARS-CoV-2",
                   "shedding=1000, symp=0.01" = "High shedding\nLow symptomatic\ne.g. Poliovirus"
)

# compute percentage where both systems detect
annot_df <- rtot %>%
  filter(profile != "shedding=0.01, symp=0.2") %>%
  group_by(profile) %>%
  summarise(
    pct_both = round(mean(detection == "both") * 100, 0)
  ) %>%
  mutate(
    label = sprintf("based on %.0f%%\nof outbreaks", pct_both)
  )

# 1. create the bar data
bars_df <- rtot %>%
  filter(profile != "shedding=0.01, symp=0.2") %>%
  group_by(profile, detection) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(profile) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 2. define position of the bars
bars_df <- bars_df %>%
  mutate(
    ymin = 0,
    ymax = prop
  ) %>%
  filter(detection != "both")

# 3. find necessary width of bars
ymax_target <- 0.02
p_max <- max(bars_df$prop)
bar_width <- p_max / ymax_target

# 4. rescale ymax and define the bar position
bars_df <- bars_df %>%
  mutate(
    height = prop / bar_width,
    xmin = case_when(
      detection == "cl" ~ -200 - bar_width,
      detection == "ww" ~ 400
    ),
    xmax = case_when(
      detection == "cl" ~ -200,
      detection == "ww" ~ 400 + bar_width
    ),
    ymax = height   # equals ymax_target
  )

# 5. calculate fraction of "both" detections per profile
both_frac <- rtot %>%
  filter(profile != "shedding=0.01, symp=0.2") %>%
  group_by(profile) %>%
  summarise(frac_both = mean(detection == "both")) %>%
  ungroup()

scale_factors <- both_frac$frac_both
names(scale_factors) <- both_frac$profile

# 6. precompute densities per profile
dens_df <- rtot %>%
  filter(detection == "both") %>%
  group_by(profile) %>%
  summarise(dens = list(density(time_saved)), .groups = "drop") %>%
  mutate(
    x = map(dens, ~ .$x),
    y = map2(dens, profile, ~ .$y * scale_factors[.y])
  ) %>%
  unnest(c(x, y))

# 7. build the plot with precomputed densities
pcv4 <- ggplot() +
  # bars
  geom_rect(
    data = bars_df,
    aes(fill = detection, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    alpha = 1,
    color = "black",
    inherit.aes = FALSE
  ) +
  # precomputed density
  geom_area(
    data = dens_df,
    aes(x = x, y = y), fill = "grey30", position = "identity"
  ) +
  scale_fill_manual(values = c( "salmon4",
                                "#5B7C99",
                                "grey20",
                                "#D3D3D3" ),
                    labels = c("Wastewater", "Clinical", "Both", "Neither"),
                    breaks = c("ww", "cl", "both", "neither"))+
  geom_line(
    data = dens_df,
    aes(x = x, y = y), color = "black"
  ) +
  facet_wrap(
    ~ profile,
    nrow = 4,
    strip.position = "right",
    labeller = labeller(profile = custom_labels)
  ) +
  geom_vline(xintercept = 0, col = "black", linetype = 2) +
  scale_x_continuous(
    breaks = seq(-100, 300, 100),
    expand = c(0, 0)
  )+
  coord_cartesian(ylim = c(0, 0.02), xlim = c(-250, 450)) +
  theme_classic() +
  theme(
    panel.spacing = unit(0, "lines"),
    axis.line = element_line(color = "black"),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(color = "white"),
    text = element_text(size = 14),
    legend.position = "none",
    axis.line.x = element_blank()
  ) +
  labs(
    x = "Days of advanced warning from wastewater\nwhen both methods detect the outbreak",
    y = "Density"
  )+
  geom_hline(yintercept = 0, color = "black")

pBCv4 <- cowplot::plot_grid(pb, pcv4, nrow = 2, rel_heights = c(0.75,1))

ggsave(pBCv4, filename = here(
  "figures/fig2BC.png"
  ), dpi = 300, width = 6, height = 9)

