# This script takes the newly reshaped Hewitt et al data and refits the --------
# detection curves

# load libraries
library(tidyverse)
library(rstan)
library(loo)

## read in clean Hewitt et al data updated to reflect fecal shedding profile
dat <- read.csv(
  here(
    my_path,
    "sc2_data_fecal_update.csv")
)
datWWTP <- dat %>% filter(!is.na(ww_detection_WWTP))

## Fit logistic regression to the data using fecal shedding profile
fit_logisticWWTP_fecal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "logistic_reg_model.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.fecal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_logisticWWTP_fecal, here(
  "results/ww_sensitivity_model_fits/refit_logistic_WWTP_fecal.rds")
  )

## Fit logistic regression to the data using nasal shedding profile
fit_logisticWWTP_nasal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "logistic_reg_model.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.nasal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_logisticWWTP_nasal, here(
  "results/ww_sensitivity_model_fits/refit_logistic_WWTP_nasal.rds")
  )

## Try fitting to log shedders as perhaps the issue is that the relationship
## between detection and shedders is not linear on the logit scale

# using the fecal shedding profile
fit_logisticWWTPlog_fecal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "logistic_reg_model_log.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.fecal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_logisticWWTPlog_fecal, here(
  "results/ww_sensitivity_model_fits/refit_logistic_WWTP_logshed_fecal.rds")
  )

## save for use in wastewatcher
tmp <- rstan::summary(fit_logisticWWTPlog_fecal, pars = c("alpha", "beta"))$summary
saveRDS(tmp, file = here(
  "results/estimates/ww_detection_params.rds")
        )

# and using the original nasal shedding profile
fit_logisticWWTPlog_nasal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "logistic_reg_model_log.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.nasal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_logisticWWTPlog_nasal, here(
  "results/ww_sensitivity_model_fits/refit_logistic_WWTP_logshed_nasal.rds")
)

## Try fitting to log shedders as perhaps the issue is that the relationship
## between detection and shedders is not linear on the logit scale

# using the fecal shedding profile
fit_hillWWTPlog_fecal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "hill_model_log.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.fecal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_hillWWTPlog_fecal, here(
  "results/ww_sensitivity_model_fits/refit_hill_WWTP_logshed_fecal.rds")
)

# and using the original nasal shedding profile
fit_hillWWTPlog_nasal <- stan(
  file = here("R", "wastewater_sensitivity_models",
              "hill_model_log.stan"),
  data = list(N = dim(datWWTP)[1],
              y = datWWTP$ww_detection_WWTP,
              N_shedders = datWWTP$effective_shedding.nasal),
  iter = 4000,
  chains = 4,
  seed = 123
)
saveRDS(fit_hillWWTPlog_nasal, here(
  "results/ww_sensitivity_model_fits/refit_hill_WWTP_logshed_nasal.rds")
)

## compare fits
fit_logistic_nasal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits/refit_logistic_WWTP_nasal.rds")
  )
fit_logistic_fecal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits/refit_logistic_WWTP_fecal.rds")
  )
fit_loglogistic_nasal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits/refit_logistic_WWTP_logshed_nasal.rds")
  )
fit_loglogistic_fecal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits/refit_logistic_WWTP_logshed_fecal.rds")
  )
fit_hill_nasal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits/refit_hill_WWTP_logshed_nasal.rds")
  )
fit_hill_fecal <- readRDS(
  here(
    "results/ww_sensitivity_model_fits//refit_hill_WWTP_logshed_fecal.rds")
  )


# 1. extract the log-likelihoods and compute LOO for each model ----------------

loo_logistic_nasal <- loo(
  extract_log_lik(fit_logistic_nasal, merge_chains = FALSE)
  )
loo_logistic_fecal <- loo(
  extract_log_lik(fit_logistic_fecal, merge_chains = FALSE)
  )
loo_loglogistic_nasal <- loo(
  extract_log_lik(fit_loglogistic_nasal, merge_chains = FALSE)
  )
loo_loglogistic_fecal <- loo(
  extract_log_lik(fit_loglogistic_fecal, merge_chains = FALSE)
  )
loo_hill_nasal <- loo(
  extract_log_lik(fit_hill_nasal, merge_chains = FALSE)
  )
loo_hill_fecal <- loo(
  extract_log_lik(fit_hill_fecal, merge_chains = FALSE)
  )

# 2. collect into list for comparison ------------------------------------------

loos <- list(
  logistic_nasal     = loo_logistic_nasal,
  logistic_fecal     = loo_logistic_fecal,
  loglogistic_nasal  = loo_loglogistic_nasal,
  loglogistic_fecal  = loo_loglogistic_fecal,
  hill_nasal         = loo_hill_nasal,
  hill_fecal         = loo_hill_fecal
)

# 3. compare loos --------------------------------------------------------------

comp <- loo_compare(loos)

# 4. save table ----------------------------------------------------------------

saveRDS(comp, file = here("results", "ww_sensitivity_model_fits", "loo_comp.rds"))

# 5. extract parameter estimates for plotting ----------------------------------

extract_pars_for_plot <- function(stanfit, N, model){

  post <- rstan::extract(stanfit)

  if(model == "hill"){
    EC50 <- exp(post$log_EC50)
    h <- post$h

    p_post <- sapply(1:1000, function(i) {
      1 / (1 + (EC50[i] / (N + epsilon))^h[i])
    })

  } else if(model == "standard"){
    p_post <- sapply(1:1000, function(i) {
      1 / (1 + exp(-(post$alpha[i] + post$beta[i] * N)))
    })

  } else if(model == "log10-transformed"){
    p_post <- sapply(1:1000, function(i) {
      1 / (1 + exp(-(post$alpha[i] + post$beta[i] * log10(N))))
    })
  }

  p_summary <- data.frame(
    N_shedders = N_seq,
    N = N,
    p_mean  = apply(p_post, 1, mean),
    p_lower = apply(p_post, 1, quantile, 0.025),
    p_upper = apply(p_post, 1, quantile, 0.975)
  )

  return(p_summary)
}

## N shedders vector for predicted probability of detection for plotting
N_seq <- seq(0, 50, length.out = 500)
epsilon <- 1e-3
N <- N_seq + epsilon

p_summary_fecal <- extract_pars_for_plot(
  stanfit = fit_logistic_fecal, N = N, model = "standard"
  )
p_summary_fecal_ll <- extract_pars_for_plot(
  stanfit = fit_loglogistic_fecal, N = N, model = "log10-transformed"
  )
p_summary_fecal_hill <- extract_pars_for_plot(
  stanfit = fit_hill_fecal, N = N, model = "hill"
  )
p_summary_nasal <- extract_pars_for_plot(
  stanfit = fit_logistic_nasal, N = N, model = "standard"
  )
p_summary_nasal_ll <- extract_pars_for_plot(
  stanfit = fit_loglogistic_nasal, N = N, model = "log10-transformed"
  )
p_summary_nasal_hill <- extract_pars_for_plot(
  stanfit = fit_hill_nasal, N = N, model = "hill"
  )

# 6. merge for plotting --------------------------------------------------------

p_summary <- bind_rows(p_summary_fecal,
                       p_summary_fecal_ll,
                       p_summary_fecal_hill,
                       p_summary_nasal,
                       p_summary_nasal_ll,
                       p_summary_nasal_hill,
                       .id = "id")
p_summary <- p_summary %>%
  mutate(
    profile = if_else(id <=3, "fecal", "nasal"),
    model = case_when((id == 1 | id == 4) ~ "standard logistic",
                      (id == 2 | id == 5) ~ "log10-transformed logistic",
                      (id == 3 | id == 6) ~ "hill"))

saveRDS(p_summary, "results/ww_sensitivity_model_fits/shedding_curves.rds")

# 7. reformat data to add as rug plot ------------------------------------------

dat_rug <- datWWTP %>%
  select(ww_detection_WWTP,
         effective_shedding.fecal,
         effective_shedding.nasal) %>%
  pivot_longer(
    cols = c(effective_shedding.fecal, effective_shedding.nasal),
    names_to = "profile",
    values_to = "N_eff"
  ) %>%
  mutate(
    profile = recode(profile,
                     "effective_shedding.fecal" = "fecal",
                     "effective_shedding.nasal" = "nasal"),
    x_rug = log10(N_eff),
    y_rug = ifelse(ww_detection_WWTP == 1, 1, 0)  # top = 1, bottom = 0
  )

# 7. plot and save -------------------------------------------------------------

## hewitt et al original curve parameters for reference
alpha_ref <- -1.21
beta_ref  <- 0.31
p_ref <- data.frame(N_shedders = N,
                    p_mean = 1 / (1 + exp(-(alpha_ref + beta_ref * N))))

## order factors
p_summary <- p_summary %>%
  mutate(model = factor(model, levels = c("log10-transformed logistic", "hill", "standard")))

## facet labels
labeller_plot <- labeller(
  model = c(
    standard = "Standard logistic",
    "log10-transformed logistic" = "Log10-transformed logistic regression",
    hill = "Hill function"
  ),
  profile = c(fecal = "Fecal shedding profile",
              nasal = "Nasal shedding profile")
)

p1 <- ggplot(p_summary%>%filter(model != "standard")) +
  # Credible interval ribbon
  geom_ribbon(aes(x = log10(N), ymin = p_lower, ymax = p_upper, fill = profile, group = profile),
              alpha = 0.3) +
  # Posterior mean curve
  geom_line(aes(x = log10(N), y = p_mean, color = profile, group = profile), size = 1.2) +
  # Distribution of number of shedders for detections
  geom_rug(
    data = dat_rug %>% filter(ww_detection_WWTP == 1),
    aes(x = x_rug, color = profile),
    sides = "t",
    alpha = 0.6
  ) +
  # Distribution of number of shedders for non-detections
  geom_rug(
    data = dat_rug %>% filter(ww_detection_WWTP == 0),
    aes(x = x_rug, color = profile),
    sides = "b",
    alpha = 0.6
  ) +
  # Reference curve overlay
  geom_line(data = p_ref, aes(x = log10(N_shedders), y = p_mean, color = "Original Hewitt Model 3"),
            size = 1.1, linetype = "dashed") +
  facet_wrap(~model, labeller = labeller_plot)+
  labs(
    x = "Effective number of shedders",
    y = "Detection probability",
    color = "Curve"
  ) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1,2),
                     labels = 10^c(-3, -2, -1, 0, 1,2))+
  scale_color_manual(values = c("brown", "turquoise4", "red"),
                     labels = c("Fecal shedding", "Nasal shedding",
                                "Original standard logistic curve"),
                     name = " ") +
  scale_fill_manual(values = c("brown", "turquoise4", "red"),
                    labels = c("Fecal shedding", "Nasal shedding",
                               "Original standard logistic curve"),
                    name = " ") +
  theme_minimal(base_size = 14)+
  theme(legend.position = "bottom")+ guides(fill = "none")
ggsave(p1, filename = here("figures/figS4.png"), dpi = 300)
