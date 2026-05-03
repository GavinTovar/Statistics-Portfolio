// Question 3 - Part 1
data {
  int<lower=1> I;                              // Number of groups
  int<lower=1> N;                              // Total number of observations
  int<lower=0> N0;                             // Number of observed y (i.e., where m = 0)
  array[N] real tilde_t;                       // Standardized time vector
  array[N0] real tilde_t_obs;                  // Observed Standardized time vector
  array[N] int<lower=1, upper=I> group;        // Group indices for each observation
  array[N0] int<lower=1, upper=I> g_obs;       // Group indices for each observed observation
  array[N0] real y_obs;                        // Observed responses
  array[N] int<lower=0, upper=1> m;            // Missingness indicator
  array[N] int<lower=-1, upper=1> treat;       // Treatment vector
  array[N0] int<lower=-1, upper=1> treat_obs;  // Observed Treatment vector
}

parameters {
  real<lower=0> sigma;   // Observation-level noise
  real betatreat;        // Coefficient on treatment
  real muA;              // Mean of alpha
  real muB;              // Mean of beta
  real<lower=0> tauA;    // Standard deviation of alpha
  real<lower=0> tauB;    // Standard deviation of beta
  real phi0;             // Logistic regression intercept for missingness
  real phi1;             // Logistic regression slope for missingness
  real phi2;             // Logistic regression coefficient for missingness
  real phi3;             // Logistic regression coefficient for missingness
  vector[I] alpha;       // Latent group-level effects
  vector[I] beta;        // Coefficient tied to time
}

model {
  // Making up priors
  tauB ~ gamma(8, 1);
  phi2 ~ normal(0.5, 1);
  phi3 ~ normal(0.5, 1);

  // Prior for alpha
  alpha ~ normal(muA, tauA);
  
  // Prior for beta
  beta ~ normal(muB, tauB);

  // Likelihood for observed y
  for (k in 1:N0) {
    y_obs[k] ~ normal(alpha[g_obs[k]] + beta[g_obs[k]] * tilde_t_obs[k] + betatreat * treat_obs[g_obs[k]], sigma);
  }

  // Likelihood for missingness
  for (k in 1:N) {
    m[k] ~ bernoulli_logit(phi0 + phi1 * alpha[group[k]] + phi2 * beta[group[k]] + phi3 * treat[group[k]]);
  }
}
