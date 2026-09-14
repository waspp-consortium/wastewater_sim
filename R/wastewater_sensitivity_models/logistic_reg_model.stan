/// Logistic regression model (as used by Hewitt et al)
// to be fit to case data and wastewater binary detection data

data {
  int<lower=0> N;               // number of observations
  int<lower=0,upper=1> y[N];    // detection outcome (0/1)
  real<lower=0> N_shedders[N];  // effective number of shedders
}

parameters {
  real alpha;   // intercept
  real beta;    // slope
}

transformed parameters {
  real epsilon = 1e-3; // small offset to handle zeros
}

model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);

  // Likelihood
  for (i in 1:N) {
    //real logN = log(N_shedders[i] + epsilon);
    y[i] ~ bernoulli_logit(alpha + beta * (N_shedders[i] + epsilon));
  }
}

generated quantities {
  vector[N] log_lik;

  // Log-likelihood for use in model fit comparisons
  for (i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] |
                        alpha + beta * (N_shedders[i] + epsilon));
  }
}
