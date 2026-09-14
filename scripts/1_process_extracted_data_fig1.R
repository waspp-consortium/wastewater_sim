# This script processes the raw extracted data from the review for use in the
# data summary in Figure 1 and the downstream modelling described in Figure 2

# it depends on:
### - data/raw/literature_review_results/Table_S1.xlsx
### - data/raw/literature_review_results/Table_symptoms_durations.csv

# it produces
### - data/processed/shedding/VL_by_study_per_unit_sample.rds
### - data/processed/shedding/VL_by_virus_per_day_ranked.rds
### - data/processed/shedding/VL_by_virus_per_infection.rds

# load libraries
library(tidyverse)
library(here)

################################################################################
# 1 - process for Figure 1 comparing VL across viruses & studies by sample type
################################################################################

# load the results of the literature review from Table S1
tab <- readxl::read_xlsx(
  here(
  "data", "raw", "literature_review_results", "Table_S1.xlsx"),
                         na = c("NA", "")
  )

# filter out papers without viral load results
# convert /wet-g and /wet-mg to /ml (Rose at al density = 1.06g/ml)
# convert dry weight --> wet weight (Rose at al 128g wet --> 29g dry)
tab <- tab %>%
  filter(
    VL_available == "Y"
    ) %>%
  mutate(
    range_central_VL_1 = as.numeric(range_central_VL_1),
    multiplier = case_when(
      VL_unit == "copies_g"  ~ 1/1.06,
      VL_unit == "copies_mg" ~ 1000/1.06,
      VL_unit == "gene_copies_mg_dry_weight" ~ (1000/(1.06 * (128/29))),
      TRUE ~ 1   # this keeps NA rows unchanged
    ),
    across(
      c(central_VL,
        range_central_VL_1, range_central_VL_2,
        variability_VL_1, variability_VL_2,
        uncertainty_mean_VL_1, uncertainty_mean_VL_2),
      ~ .x * multiplier
    )
  ) %>%
  select(
    -multiplier
    )


# re-order to maintain virus family structure
virus_order <- c("LASV", "MERS-CoV", "SARS-CoV-2", "EBOV",
                 "DENV", "YFV", "ZIKV", "CCHFV", "IAV", "SFTSV",
                 "MPXV", "CHIKV", "PV")

# categorise by = or > 1 individual contributing to measure
tab <- tab%>%
  mutate(
    virus = factor(
      virus, levels = virus_order
      ),
    indiv = case_when(
      n_central_VL==1 ~ "single",
      n_central_VL>1 ~ "multi",
      # unclear how many people Bozza et al tested, but >1
      virus == "ZIKV" & study == "Bozza et al 2019" ~ "multi",
      is.na(n_central_VL) ~ as.character(NA)))

# aggregate Thielebein et al timepoints --> range across time
thieb <- tab %>%
  filter(
    study == "Thielebein et al 2022"
    )%>%
  group_by(
    study, virus, family, sample_type
    )%>%
  summarise(
    variability_VL_1 = min(central_VL),
    variability_VL_2 = max(central_VL),
    variability_type = "range",
    variability_across = "time",
    n_central_VL = 110,
    individuals_type = "confirmed_convalescent"
    )

tab <- bind_rows(
  tab %>%
    filter(
      study != "Thielebein et al 2022"
      ),
  thieb
  )

# avoid double counting:
# remove To et al as covered by Lowry et al meta-analysis
tab <- tab %>%
  filter(!(study == "To et al 2009" & sample_type == "Stool"))

saveRDS(tab,
        file = here(
          "data/processed/shedding/VL_by_study_per_unit_sample.rds"
          )
        )
################################################################################
# 2 - process study specific VL/unit sample --> virus specific VL/day
################################################################################

tab <- readRDS(
  here(
    "data/processed/shedding/VL_by_study_per_unit_sample.rds"
    )
  )

# first get to a summary measure for urine and stool for each virus
## filter out Wannigama JN.1 and BA.2.86 as these shed much more and were
## not present when Hewitt et al conducted their sensitivity estimates.
tab1 <- tab %>%
  filter( # filter out studies not reporting central/variability estimate
    !is.na(central_VL) |
      !is.na(variability_VL_1) |
      !is.na(uncertainty_mean_VL_1),
    study != "Wannigama et al 2024 including BA.2.86 & JN.1"
    ) # filter out extra Wannigama value with later variants to avoid duplication

# compute sample-size-weighted geometric means for each panel
## weight look-up table
weights <- tab1 %>%
  group_by(
    virus, sample_type
  ) %>%
  summarise(
    tot_N = sum(n_central_VL)
  )
## pull weights
tab2 <- left_join(
  tab1, weights,
  by = c(
    "virus", "sample_type"
  )
)%>%
  mutate( # get values to take geometric mean from (log scale)
    x_value = case_when(
      !is.na(central_VL) ~ log10(central_VL),
      !is.na(variability_VL_1) & !is.na(variability_VL_2) ~
        (log10(variability_VL_1) + log10(variability_VL_2)) / 2,
      TRUE ~ NA_real_
    ),
    weight = n_central_VL/tot_N
  )%>%
  group_by(
  sample_type, virus
)%>%
  summarise(
    overall_lowbound = if (first(virus) == "SARS-CoV-2") {
      min(central_VL, na.rm = TRUE)
    } else {
      min(
        c(central_VL, variability_VL_1, uncertainty_mean_VL_1),
        na.rm = TRUE
      )
    },

    overall_highbound = if (first(virus) == "SARS-CoV-2") {
      max(central_VL, na.rm = TRUE)
    } else {
      max(
        c(central_VL, variability_VL_2, uncertainty_mean_VL_2),
        na.rm = TRUE
      )
    },

    overall_central = if (first(virus) == "ZIKV") {
      10^(mean(x_value, na.rm = TRUE))
    } else {
      10^(sum(x_value * weight, na.rm = TRUE))
    },

    n_indiv = sum(n_central_VL, na.rm = TRUE),
    n_studies = n(),
    .groups = "drop"
  )

# since we do not have an estimate of variability or uncertainty for MERS-CoV
# we will present a conservative range of 100x more/less than the central value
tab2 <- tab2%>%
  mutate(
    overall_lowbound = ifelse(virus == "MERS-CoV",
      overall_central/100,
      overall_lowbound
    ),
    overall_highbound = ifelse(virus == "MERS-CoV",
      overall_central*100,
      overall_highbound
    )
  )

# since we do not have an estimate of variability or uncertainty for LASV in stool
# we will present a conservative range of 100x more/less than the central value
tab2 <- tab2%>%
  mutate(
    overall_lowbound = ifelse(virus == "LASV" & sample_type == "Stool",
                              overall_central/100,
                              overall_lowbound
    ),
    overall_highbound = ifelse(virus == "LASV" & sample_type == "Stool",
                               overall_central*100,
                               overall_highbound
    )
  )

# we do not have a measure of viral copy numbers of MPXV in stool. However we do
# have 3 studies that detected MPXV in stool (often more commonly than in urine),
# and measured Ct values for both stool and urine.
# The first reports Ct ~30 in stool and Ct<30 in urine (Yang et al 2024b),
# the second has only presented Ct values in a figure and it is hard to extract
# them (Guo et al 2024), and the third presents a table of Ct values (n=12):
# https://doi.org/10.2807/1560-7917.ES.2022.27.28.2200503

stoolCt_mpxv_peiro_mestres <- c(23.4, 24.0, 21.5, 19.9, 24.7, 26.4, 28.3, 31.4,
                              30.2, 29.8, 20.0, 23.4, 20.9, 17.8)
urineCt_mpxv_peiro_mestres <- c(37.3, 27.0, 39.1, 39.3, 39.4, 24.4, 19.1, 40,
                                40, 31.9, 35.0, 29.0, 26.7)
diffCt <- mean(urineCt_mpxv_peiro_mestres) - mean(stoolCt_mpxv_peiro_mestres)

## viral shedding difference is therefore ~ 2^diffCt

## so we shall assume viral load in stool of 370 * urine +/- 100*
mpxv_stool <- data.frame(virus = "MPXV", sample_type = "Stool",
                         overall_central = (2^diffCt) * tab2 %>%
                           filter(virus == "MPXV", sample_type == "Urine")%>%
                           pull(overall_central),
                         n_indiv = 0, n_studies = 0)%>%
  mutate(overall_lowbound = overall_central/100,
         overall_highbound = overall_central*100)
tab2 <- bind_rows(tab2, mpxv_stool)

# The percentage of people shedding is tiny (ref) so we will assume no contribution
# of shedding in urine for IAV
# tab1 <- tab1 %>% filter(!(virus == "IAV" & sample_type == "Urine"))


# convert shedding estimates from /sample type /ml to shedding /day

## specify parameters -----------------------------------------------------------
stool_density <- 1.06 # grams/mL
## ref: Penn et al 2018 https://doi.org/10.1016/j.watres.2017.12.063
stool_vol <- c("low"=72, "upp"=470, "median"=128)/stool_density # mL/person/day
## refs: Cummings et al 1992 https://doi.org/10.1016/0016-5085(92)91435-7
##       Rose et al 2015 https://doi.org/10.1080/10643389.2014.1000761
urine_vol <- c("low"=600, "upp"=2600, "median"=1420) # mL/person/day
## ref: Rose et al 2015 https://doi.org/10.1080/10643389.2014.1000761


# Multiply by daily waste volumes ----------------------------------------------
tab3 <- tab2 %>%
  mutate(overall_lowbound = if_else(sample_type == "Stool",
                                    overall_lowbound*(stool_vol["low"]),
                                    overall_lowbound*urine_vol["low"]),
         overall_highbound = if_else(sample_type == "Stool",
                                     overall_highbound*(stool_vol["upp"]),
                                     overall_highbound*urine_vol["upp"]),
         overall_central = if_else(sample_type == "Stool",
                                   overall_central*(stool_vol["median"]),
                                   overall_central*urine_vol["median"]))
# group by virus, sum over stool & urine and save for ranking by daily shedding
tab4 <- tab3 %>%
  group_by(virus)%>%
  summarise(vl_low = sum(overall_lowbound),
            vl_upp = sum(overall_highbound),
            vl_average = sum(overall_central)
  )%>%
  arrange(desc(vl_average))
saveRDS(tab4, here("data/processed/shedding/VL_by_virus_per_day_ranked.rds"))


################################################################################
# 3 - combine with duration --> VL shed per infection relative to SARS-CoV-2
################################################################################

# merge with duration and percentage symptomatic for quick figure
pars_mod <- read.csv(here("data/raw/literature_review_results/duration_symptomprop.csv"))
dat <- left_join(tab3, pars_mod, by = c("virus", "sample_type"))

# 1. There is not enough data to make any rational assumptions about nipah virus
#    or hantaviruses or pheniuviruses or MERS-CoV so we will remove these.
#    We will also remove YFV as shedding data is only available for the vaccine/
#    convalescents
dat <- dat %>% filter(!virus %in% c("NIV", "HTNV", "SFTSV", "MERS-CoV", "YFV"))

# 3. We do not have a duration of shedding in stool for Lassa virus but we know
#    it is present (day 11, mean 10,000 copies/mL which is 10x mean in urine)
# we will assume that the range of shedding duration in stool is 1 week as a
# conservative lower bound, and the upper is the same as the upper for urine
dat[which(dat$virus == "LASV" & dat$sample_type == "Stool"), "duration_l"] <- 1
dat[which(dat$virus == "LASV" & dat$sample_type == "Stool"), "duration_u"] <-
  dat[which(dat$virus == "LASV" & dat$sample_type == "Urine"), "duration_u"]

# 4. multiply daily VL by duration in days
dat <- dat %>%
  mutate(vl_low = case_when(!is.na(duration_l) ~ overall_lowbound*duration_l*7,
                            is.na(duration_l) ~ overall_lowbound*average.duration*7,
                            TRUE ~ NA),
         vl_upp = case_when(!is.na(duration_u) ~ overall_highbound*duration_u*7,
                            is.na(duration_u) ~ overall_highbound*average.duration*7,
                            TRUE ~ NA),
         vl_average = overall_central*average.duration*7)

# 5. group by virus
df <- dat %>%
  group_by(virus, perc_symp, perc_symp_l, perc_symp_u)%>%
  summarise(vl_low = sum(vl_low),
            vl_upp = sum(vl_upp),
            vl_average = sum(vl_average)
  )
# 6. add family
df$family <- c("Nairoviruses", "Togaviruses", "Flaviviruses",
                 "Filoviruses", "Orthomyxoviruses", "Arenaviruses",
                  "Poxviruses", "Picornaviruses",
                 "Coronaviruses", "Flaviviruses")
# 7. save
saveRDS(df, file = here("data/processed/shedding/VL_by_virus_per_infection.rds"))
write.csv(df, file = here("data/processed/shedding/VL_by_virus_per_infection.csv"))
