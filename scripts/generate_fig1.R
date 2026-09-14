# This script plots the processed VL data by study per unit sample (Figure 1)

# it depends on:
### - data/processed/shedding/VL_by_study_per_unit_sample.rds - from script 1
### - data/processed/shedding/VL_by_virus_per_day_ranked.rds - from script 1

# it produces
### - figures/fig1.png

# load libraries
library(tidyverse)
library(tidytext) # for reordering factors neatly to aid plotting y axis

################################################################################
#
################################################################################

# load processed shedding data
tab <- readRDS(
  here(
    "data", "processed", "shedding", "VL_by_study_per_unit_sample.rds"
    )
  )

# get rid of those that only reported a maximum viral load
tab2 <- tab %>%
  filter(!is.na(central_VL) | !is.na(range_central_VL_1) |
           !is.na(variability_VL_1))

# combine meta-analyses variability and uncertainty for joint plotting
tab2 <- tab2 %>%
  mutate(
    lower_bound = case_when(
      meta_analysis == "Y" & !is.na(
        uncertainty_mean_VL_1
        ) ~ uncertainty_mean_VL_1,
      meta_analysis == "Y" & is.na(
        uncertainty_mean_VL_1
        ) ~ variability_VL_1,
      is.na(
        meta_analysis
        ) ~ NA,
      TRUE ~ NA
    ),
    upper_bound = case_when(
      meta_analysis == "Y" & !is.na(
        uncertainty_mean_VL_2
        ) ~ uncertainty_mean_VL_2,
      meta_analysis == "Y" & is.na(
        uncertainty_mean_VL_2
        ) ~ variability_VL_2,
      is.na(
        meta_analysis
        ) ~ NA,
      TRUE ~ NA))

# create variable for measure saying what bounds describe
# (variability/uncertainty etc)
tab2 <- tab2 %>%
  mutate(
    measure = case_when(
      !is.na(variability_type) ~ paste0(
        variability_type, " across\n", variability_across
        ),
      is.na(
        variability_type
        ) & (
          study == "Overall"
          ) ~ "range across\nstudies",
      is.na(
        variability_type
        ) & !is.na(
          range_central_VL_1
          ) ~ "range across\ntime PSO",
      TRUE~as.character(
        NA
        )
      )
    )

# create pretty names and include sample size, n
tab3 <- tab2 %>%
  mutate(
    study_clean = str_remove(study, "___.*$")
    ) %>%
  group_by(
    study_clean
    ) %>%
  mutate(
    n_u = n_central_VL[sample_type == "Urine"][1],
    n_s = n_central_VL[sample_type == "Stool"][1]
  ) %>%
  ungroup() %>%
  mutate(
    # build each component only if it exists
    part_u = ifelse(!is.na(n_u),
                    paste0("n[u]==", n_u),
                    NA),

    part_s = ifelse(!is.na(n_s),
                    paste0("n[s]==", n_s),
                    NA),
    # combine only non-NA parts
    n_part = apply(
      cbind(part_u, part_s),
      1,
      function(x) {
        x <- na.omit(x)
        if (length(x) == 0) return(NA)
        paste(x, collapse = "*','~")
      }
    ),
    # final label
    study_label = ifelse(
      is.na(n_part),
      paste0("'", study_clean, "'"),
      paste0(
        "'", study_clean, "'",
        "~'['*", n_part, "*']'"
      )
    )
  )

# add asterisk to Taniuchi et al as we only know the number of samples, not
# how many infants they came from (up to 270)
tab3 <- tab3 %>%
  mutate(
    study_label = ifelse(
      study_clean == "Taniuchi et al 2014",
      "'Taniuchi et al 2014'~'['*n[s]==237*'*'*']'",
      study_label
    )
  )

# manual edit of Bozza for ZIKV (it has otherwise inherited the n for CHIKV)
tab3 <- tab3 %>%
  mutate(
    study_label = ifelse(
      study_clean == "Bozza et al 2019" & virus == "ZIKV",
      paste0(
        "'Bozza et al 2019'~'['*n[u]=='NA'*'**'*']'"
      ),
      study_label
    )
  )
# manual edit for Wannigama excluding BA.2.86 & JN.1 (over two lines nicely)
tab3 <- tab3 %>%
  mutate(
    study_label = ifelse(
      study_clean == "Wannigama et al 2024 excl. BA.2.86 & JN.1",
      "atop('Wannigama et al 2024 excl.', 'BA.2.86 & JN.1 '~'['*n[s]=='72'*']')",
      study_label
    )
  )
# compute sample-size-weighted geometric means for each panel
## weight look-up table
weights <- tab3 %>%
  filter(
    study != "Wannigama et al 2024 including BA.2.86 & JN.1"
  )%>%
  group_by(
    virus, sample_type
    ) %>%
  summarise(
    tot_N = sum(n_central_VL)
    )
## pull weights
tabx <- left_join(
  tab3, weights,
  by = c(
    "virus", "sample_type"
    )
  )
## compute geometric means
panel_means <- tabx %>%
  filter(
    study != "Wannigama et al 2024 including BA.2.86 & JN.1"
  )%>%
  mutate(
    x_value = case_when(
      !is.na(central_VL) ~ log10(central_VL),
      !is.na(variability_VL_1) & !is.na(variability_VL_2) ~
        (log10(variability_VL_1) + log10(variability_VL_2)) / 2,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    weight = n_central_VL/tot_N
  )%>%
  group_by(
    virus, sample_type
    ) %>%
  summarise(
    panel_mean = if_else(
      first(virus) == "ZIKV",
      mean(x_value, na.rm = T),
      sum(x_value*weight, na.rm = T)),
    .groups = "drop"
  )

# place blank rectangles over the empty panels
empty_panels <- expand_grid(
  virus = unique(
    tab3$virus
    ),
  sample_type = unique(
    tab3$sample_type
    )
) %>%
  anti_join(
    tab3, by = c("virus","sample_type")
    )

# we want to order by overall daily magnitude shed
rank <- readRDS(
  here(
    "data","processed","shedding","VL_by_virus_per_day_ranked.rds"
    )
)
virus_order <- rank$virus

tab3 <- tab3 %>%
  mutate(
    virus = factor(
      virus, levels = virus_order
      ),
    individuals_type = factor(
      individuals_type, levels = c("confirmed_case",
                                   "confirmed_convalescent",
                                   "vaccinated",
                                   "fyi")
    )
    )

# plot
p3 <- ggplot(
  tab3
  )+
  geom_rect(
    data = empty_panels,
    aes(
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ),
    inherit.aes = FALSE,
    fill = "grey95"
  )+
  geom_errorbar(
    data = tab3 %>% filter(
      variability_across!="studies"),
                aes(
                  y = study_label,
                  xmin = log10(variability_VL_1),
                  xmax = log10(variability_VL_2),
                  col = individuals_type
                  ),
    alpha = 0.75, linewidth = 3, width = 0,
    position = position_dodge(width = 0.5)
    )+
  geom_errorbar(
    aes(
      y = study_label,
      xmin = log10(lower_bound),
      xmax = log10(upper_bound),
      col = individuals_type),
    width = 0.2,
    linewidth = 1,
    position = position_dodge(width = 0.5)
    )+
  geom_point(
    aes(
      y = study_label,
      x = log10(central_VL),
      shape = indiv,
      fill = individuals_type),
    position = position_dodge(width = 0.5),
    size = 3)+
  facet_grid(
    rows = vars(virus),
    cols = vars(sample_type),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
               virus = c("MERS-CoV" = "MERS\nCoV",
                         "SARS-CoV-2" = "SARS\nCoV-2")
               )
    )+
  scale_colour_manual(values = c(  "#4E6A8C",
                                   "#C46A4A",
                                   "#5A8F6B",
                                   "pink2"),
                      name = "",
                      labels = c("Acute cases",
                                 "Convalescent cases",
                                 "Vaccinated shedders",
                                 "Excluded from mean"),
                      na.value = "#8A7D6A")+
  scale_fill_manual(values = c(  "#4E6A8C",
                                 "#C46A4A",
                                 "#5A8F6B",
                                 "pink2"),
                    name = "", guide = "none",
                    na.value = "#8A7D6A")+
  scale_shape_manual(values = c(24,21),
                     labels = c(
                       "average across\nn individuals",
                       "average for\nsingle individual"
                       ),
                     name = "",
                     na.translate = FALSE)+
  xlab(
    "Viral load (Log10 copies/ml)"
    )+
  ylab(
    ""
    )+
  guides(
    colour = guide_legend(
      order = 1,
      ncol = 2
    ),
    shape = guide_legend(
      order = 2
    )
  )+
  geom_vline(
    data = panel_means,
    aes(xintercept = panel_mean),
    linetype = "dotted",
    linewidth = 0.8
  )+
  theme_bw()+
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.box = "vertical",
    legend.margin = margin(-5, 0, 0, 0),

    axis.text.y = element_text(size = 12, colour = "black",  lineheight = 0.7),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),

    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    strip.text.y = element_text(size = 8),

    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )+
  scale_y_reordered(
    labels = function(x) parse(text = x)
    )

ggsave(p3, filename = "figures/fig1.png",
       dpi = 300, width = 11, height = 11, unit = "in")

