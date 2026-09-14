# This script generates Figure 2A

# It depends on:
# - results of "5_simulations_fig2.R"

# It produces:
# - "figures/fig2A.png"

# Load libraries
library(tidyverse)
library(cowplot)
library(purrr)
library(ggnewscale)
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

## reshape - group by parameter set
ptot <- rtot %>%
  dplyr::mutate(
    detected = if_else(
      !is.na(day1WW) | !is.na(day1Clinic), "Y", "N"
      ),
    first_detect_WW = case_when(
      (!is.na(day1WW) & !is.na(day1Clinic) & day1WW<day1Clinic) |
        (!is.na(day1WW) & is.na(day1Clinic)) ~ "Y",
      (is.na(day1WW) & !is.na(day1Clinic)) |
        ((!is.na(day1WW) & !is.na(day1Clinic)) & day1Clinic<=day1WW) ~  "N",
      (is.na(day1WW) & is.na(day1Clinic)) ~ as.character(NA),
      TRUE ~ as.character(NA)
      )
    )%>%
  group_by(
    symptoms, prob_diagnosis, delay_diagnosis_shape, shed_rel_sc2
    )%>%
  dplyr::summarise(
    n = n(),
    percentage_detected_WW = 100 * sum(
      !is.na(day1WW)
      )/n,
    percentage_detected_clinical = 100 * sum(
      !is.na(day1Clinic)
      )/n,
    percentage_WW_best = 100 * sum(
      first_detect_WW == "Y", na.rm = T
      )/sum(
        detected == "Y", na.rm= T
        ),
    median_delay_WW = median(day1WW, na.rm = T),
    median_delay_clinical = median(day1Clinic, na.rm = T),
    percentage_first_detect_WW = 100*sum(
      first_detect_WW == "Y", na.rm = T)/sum(
        detected == "Y"
        ),
    majority_first_detect_WW = if_else(
      sum(
        first_detect_WW == "Y", na.rm = T
        ) > sum(
          first_detect_WW == "N", na.rm = TRUE
          ), "Y", "N")
    )%>%
  dplyr::mutate(
    WW_more_likely_detect = if_else(
      percentage_detected_WW>percentage_detected_clinical, "Y", "N"
      )
    )

## Look at where WW adds valuer (EITHER by detecting more often OR
## by detecting earlier in cases where both WW and clinical detect).
## (Also looked at where the line is where WW is best in 95% of
## detected outbreaks and where WW is best in only 5% of the detected outbreaks)

  thresholds <- ptot %>%
    dplyr::group_by(shed_rel_sc2, prob_diagnosis, delay_diagnosis_shape) %>%
    dplyr::summarise(
      threshold_symptoms = max(
        ifelse(majority_first_detect_WW == "Y", symptoms, 0),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  thresholds05 <- ptot %>%
    group_by(shed_rel_sc2, prob_diagnosis, delay_diagnosis_shape) %>%
    summarise(
      threshold_symptoms = max(
        ifelse(percentage_first_detect_WW >=5, symptoms, 0),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  thresholds95 <- ptot %>%
    group_by(shed_rel_sc2, prob_diagnosis, delay_diagnosis_shape) %>%
    summarise(
      threshold_symptoms = max(
        ifelse(percentage_first_detect_WW >=95, symptoms, 0),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  thresholds_merge <- bind_rows(
    list(thresholds05, thresholds95),
    .id = "cutoff"
  ) %>%
    mutate(
      cutoff = recode(
        cutoff,
        "1" = "5",
        "2" = "95"
      )
    )
## to extend to cover polio visually
dummy <- data.frame(shed_rel_sc2 = rep(1e+05, 3),
                    prob_diagnosis = c(0.01, 0.1, 1),
                    delay_diagnosis_shape = 6,
                    threshold_symptoms = 1)
thresholds <- bind_rows(dummy, thresholds)

## add known viruses

## read in data
vlDat <- readRDS(here(
  "data/processed/shedding/VL_by_virus_per_infection.rds")
  )
## manually add border type based on data quality, where meta means it was the
## result of a meta-analysis, var means it is the result of a single study
## which reported a measure of variability, and none means no variability was
## reported and estimates are the result of a point estimate with conservative
## variability 100x lower and higher than this is presented to avoid presenting
## faux-certainty data
vlDat$low_upp_type <- c("var", "var", "var", "none", "meta",
                        "none", "none", "var", "meta", "var")

# plot

# Base plot + shading + lines
base <- ggplot() +
  scale_x_log10(
    limits = c(1e-6, 1e+5),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_ribbon(data = thresholds %>% filter(delay_diagnosis_shape == 6),
              aes(    x = shed_rel_sc2,
                      y = 100*threshold_symptoms,
                      group = as.factor(prob_diagnosis),
                      ymin = 0, ymax = 100*threshold_symptoms),
              alpha = 0.5, col = NULL, fill = "grey80"
  ) +
  geom_rect(
    data = vlDat,
    aes(
      xmin = vl_low / as.numeric(vlDat[vlDat$virus == "SARS-CoV-2", "vl_average"]),
      xmax = vl_upp / as.numeric(vlDat[vlDat$virus == "SARS-CoV-2", "vl_average"]),
      ymax = perc_symp_l,
      ymin = perc_symp_u,
      fill = family,
      colour = family,
      linetype = low_upp_type
    ),
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "#2F3E75",  # deep indigo
      "#3B8EA5",  # muted teal
      "#4C6A2F",  # olive green
      "#C2A878",
      "#7A5195",  # dusty violet
      "#C17C3A",  # soft burnt orange
      "#6B8E23",  # moss green
      "#8C3B3B",  # muted brick
      "#5E6472"   # slate grey
    ),
    name = "Family"
  ) +
  scale_color_manual(
    values = c(
      "#2F3E75",  # deep indigo
      "#3B8EA5",  # muted teal
      "#4C6A2F",  # olive green
      "#C2A878",
      "#7A5195",  # dusty violet
      "#C17C3A",  # soft burnt orange
      "#6B8E23",  # moss green
      "#8C3B3B",  # muted brick
      "#5E6472"   # slate grey
    ),
    name = "Family"
  ) +
  scale_linetype_manual(
    values = c(
      "none" = "blank",
      "meta"  = "solid",
      "var"  = "dashed"
    ),
    name = "Bounds",
    guide= "none"
  )+ # add lines ontop of threshold shading
  ggnewscale::new_scale("linetype") +   # reset line type scale
  geom_line(data = thresholds %>% filter(delay_diagnosis_shape == 6),
            aes(
              x = shed_rel_sc2,
              y = 100 * threshold_symptoms,
              linetype = as.factor(prob_diagnosis)
            ),
            lwd = 1
  ) +
  scale_linetype_manual(
    name = "Probability of diagnosis given clinical symptoms",
    values = c(2, 1, 3)
  )+
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 14)
  )+
  xlab("Mean viral load shed (copies/day) relative to SARS-CoV-2")+
  ylab("Percentage of infections with clinical symptoms (%)")


p <- base +
  guides(
    fill   = guide_legend(nrow = 3, byrow = TRUE),
    colour = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(keywidth = 2)
  )

# Extract legends
legend_top <- get_legend(p + guides(linetype = "none"))
legend_bottom <- get_legend(p + guides(fill = "none", colour = "none"))

# Remove legends from main plot
p_clean <- p + theme(legend.position = "none")

# Assemble legends in a neat way
final_plot <- plot_grid(
  legend_top,
  p_clean,
  legend_bottom,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.12)
)

ggsave(final_plot, filename = here("figures/fig2A.png"), dpi = 300,
       width = 7, height = 9, unit = "in")
