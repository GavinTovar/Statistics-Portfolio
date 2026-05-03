// Question 3 - Part 2
data {
  int<lower=1> I;                          // Number of groups
  int<lower=1> N;                          // Total number of observations
  int<lower=0> N0;                         // Number of observed y (i.e., where m = 0)
  array[N] real tilde_t;                   // Standardized time vector
  array[N] int<lower=1, upper=I> group;    // Group indices for each observation
  array[N0] int<lower=1, upper=I> g_obs;   // Group indices for each observed observation
  array[N0] real y_obs;                    // Observed responses
  array[N] int<lower=0, upper=1> m;        // Missingness indicator
  array[N] int<lower=-1, upper=1> treat;   // Treatment vector
}

parameters {
  vector[I] alpha;       // Latent group-level effects
  vector[I] beta;        // Coefficient tied to time
  real phi2;             // Logistic regression slope coeffcient 2 for missingness
  real phi3;             // Logistic regression slope coeffcient 3 for missingness
}

model {
  // Set Hyperparameters
  real betatreat = -0.5;        // Coefficient on treatment
  real sigma = 0.25;            // Observation-level noise
  real tauA = 1;                // Standard deviation of alpha
  real tauB = sqrt(8);          // Standard deviation of beta
  real muA = 4;                 // Mean of alpha
  real muB = -5;                // Mean of beta
  real phi0 = -2;               // Logistic regression intercept for missingness
  real phi1 = 0.25;             // Logistic regression slope for missingness
  
  // Priors
  // sigma ~ student_t(3, 0, 5);
  // phi0 ~ student_t(3, 0, 1);
  // phi1 ~ student_t(3, 0, 1);
  
  // Chosen Priors
  // tauA ~ student_t(3, 0, 2);
  // tauB ~ student_t(3, 0, 2);
  phi2 ~ student_t(3, 0, 1);
  phi3 ~ student_t(3, 0, 1);
  
  // Prior for alpha
  alpha ~ normal(muA, tauA);
  
  // Prior for beta
  beta ~ normal(muB, tauB);

  // Likelihood for observed y
  for (k in 1:N0) {
    y_obs[k] ~ normal(alpha[g_obs[k]] + beta[g_obs[k]] * tilde_t[k] + betatreat * treat[g_obs[k]], sigma);
  }

  // Likelihood for missingness
  for (k in 1:N) {
    m[k] ~ bernoulli_logit(phi0 + phi1 * alpha[group[k]] + phi2 * beta[group[k]] + phi3 * treat[group[k]]);
  }
}
