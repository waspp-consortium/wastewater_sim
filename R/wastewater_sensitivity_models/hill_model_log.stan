// Hill function (log stabilised) to be fit to
// case data and wastewater binary detection data

data {
  int<lower=0> N;               // number of observations
  int<lower=0,upper=1> y[N];    // detection outcome (0/1)
  real<lower=0> N_shedders[N];  // effective number of shedders
}

parameters {
  real log_EC50;    // log half-max number of shedders
  real<lower=0> h;  // Hill coefficient (steepness of slope)
}

transformed parameters {
  real epsilon = 1e-3; // small offset to handle zeros
  real EC50 = exp(log_EC50);
}

model {
  // Priors
  log_EC50 ~ normal(0, 5);  // wide prior on the log-scale
  h ~ normal(1, 1);

  // Likelihood
  for (i in 1:N) {
    real p = 1 / (1 + pow(EC50 / N_shedders[i] + epsilon, h));
    y[i] ~ bernoulli(p);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    real p = 1 / (1 + pow(EC50 / N_shedders[i] + epsilon, h));
    log_lik[i] = bernoulli_lpmf(y[i] | p);
  }
}
