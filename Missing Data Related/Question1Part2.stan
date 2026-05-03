// Question 1 - Part 2
data {
  int<lower=1> I;                         // Number of groups
  int<lower=1> N;                         // Total number of observations
  int<lower=0> N0;                        // Number of observed y (i.e., where m = 0)
  array[N] int<lower=1, upper=I> group;   // Group indices for each observation
  array[N0] int<lower=1, upper=I> g_obs;  // Group indices for each observed observation
  array[N0] real y_obs;                   // Observed responses
  array[N] int<lower=0, upper=1> m;       // Missingness indicator
}

parameters {
  real mu;               // Hyperparameter for alpha
  real<lower=0> tau;     // Standard deviation of alpha
  real<lower=0> sigma;   // Observation-level noise
  real phi0;             // Logistic regression intercept for missingness
  real phi1;             // Logistic regression slope for missingness
  vector[I] alpha;       // Latent group-level effects
}

model {
  // Fix prior for alpha
  alpha ~ normal(mu, tau);

  // Likelihood for observed y
  for (k in 1:N0) {
    y_obs[k] ~ normal(alpha[g_obs[k]], sigma);
  }

  // Likelihood for missingness
  for (k in 1:N) {
    m[k] ~ bernoulli_logit(phi0 + phi1 * alpha[group[k]]);
  }
}
