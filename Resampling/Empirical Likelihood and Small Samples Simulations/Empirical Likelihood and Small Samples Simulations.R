### Empirical Likelihood and Small Samples Simulations
# Preliminary Syntax Check
parse("./Resampling Code/Empirical Likelihood and Small Samples Simulations.R")

library(Rcpp)
# C++ Functions
sourceCpp("./Resampling Code/Small_Sims_Code.cpp")

## Distributions
# Exponential
# Normal
# Lognormal
# Mixture of Normals
# Student's t
# Triangular
# Uniform
distributions <- c("Exponential", "Gaussian", "Lognormal",
                   "Mixture", "Student's t", "Triangular", "Uniform")

## Methods 
# z-interval                            (in C++ script)
# t-interval                            (in C++ script)
# Percentile bootstrap                  (in C++ script)
# Bias-corrected accelerated bootstrap  (in C++ script)
# Bootstrap t                           (in C++ script)
# Short bootstrap t                     (in C++ script)
# Short percentile bootstrap            (in C++ script)
# Empirical Likelihood (Chi-square)     (in C++ script)
# Empirical Likelihood (F)              (in C++ script)
methods <- c("z-interval", "t-interval", "Percentile",
             "BCa", "t", "Short t", "Short Percentile",
             "Empirical (Chi)", "Empirical (F)")

## Simulation Settings
alpha <- 0.05
alphaL <- 400
nvec <- seq(3, 20, by = 1)
# nvec <- seq(3, 8, by = 1)
# nvec <- seq(9, 14, by = 1)
# nvec <- seq(15, 20, by = 1)
nsim <- 50
bsim <- 1000
start <- proc.time() # Time simulations
# nsim x alpha grid, lower/upper/captured/length/alpha x distribution x method x sample size
int_array <- array(NA, dim = c(nsim, alphaL, 5, 7, 9, length(nvec)))

# Coverage Array Fabrication
Coverage_tables <- array(NA, dim = c(nsim, length(methods), length(methods), length(distributions), length(nvec)))
# Indexing n
N <- 1
for (n in nvec) {
  for (i in 1:nsim) {
    ### Generate data
    # Exponential
    Exponential <- rexp(n, rate = 1)
    
    # Normal
    Gaussian <- rnorm(n, mean = 0, sd = 1)
    
    # Lognormal
    if (i == 1) {
      sdlog <- sqrt(log(4.67 / 1.65^2 + 1))
      meanlog <- log(1.65) - 0.5 * sdlog^2
    }
    Lognormal <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
    
    # Mixture of Normals
    Ind <- rbinom(n, 1, p = 0.25)
    Mixture <- Ind * rnorm(n, mean = -3, sd = 1) + (1 - Ind) * rnorm(n, mean = 1, sd = 1)
    
    # Student's t
    Student <- rt(n, df = 4)
    
    # Triangular
    Tri <- runif(n, min = 0, max = 1)
    Triangular <- 2 * Tri - (Tri)^2
    
    # Uniform
    Uniform <- runif(n, min = 0, max = 1)
    
    ### Construct confidence intervals
    # z-intervals
    int_array[i, , , 1, 1, N] <- z_interval_grid(Exponential, parm = 1, alpha_length = alphaL)
    int_array[i, , , 2, 1, N] <- z_interval_grid(Gaussian, parm = 0, alpha_length = alphaL)
    int_array[i, , , 3, 1, N] <- z_interval_grid(Lognormal, parm = 1.63, alpha_length = alphaL)
    int_array[i, , , 4, 1, N] <- z_interval_grid(Mixture, parm = 0, alpha_length = alphaL)
    int_array[i, , , 5, 1, N] <- z_interval_grid(Student, parm = 0, alpha_length = alphaL)
    int_array[i, , , 6, 1, N] <- z_interval_grid(Triangular, parm = 2/3, alpha_length = alphaL)
    int_array[i, , , 7, 1, N] <- z_interval_grid(Uniform, parm = 1/2, alpha_length = alphaL)
    
    # t-intervals
    int_array[i, , , 1, 2, N] <- t_interval_grid(Exponential, parm = 1, alpha_length = alphaL)
    int_array[i, , , 2, 2, N] <- t_interval_grid(Gaussian, parm = 0, alpha_length = alphaL)
    int_array[i, , , 3, 2, N] <- t_interval_grid(Lognormal, parm = 1.63, alpha_length = alphaL)
    int_array[i, , , 4, 2, N] <- t_interval_grid(Mixture, parm = 0, alpha_length = alphaL)
    int_array[i, , , 5, 2, N] <- t_interval_grid(Student, parm = 0, alpha_length = alphaL)
    int_array[i, , , 6, 2, N] <- t_interval_grid(Triangular, parm = 2/3, alpha_length = alphaL)
    int_array[i, , , 7, 2, N] <- t_interval_grid(Uniform, parm = 1/2, alpha_length = alphaL)
    
    # Percentile bootstrap
    int_array[i, , , 1, 3, N] <- perc_boot_interval_grid(Exponential, parm = 1, B = bsim, alpha_length = alphaL)
    int_array[i, , , 2, 3, N] <- perc_boot_interval_grid(Gaussian, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 3, 3, N] <- perc_boot_interval_grid(Lognormal, parm = 1.63, B = bsim, alpha_length = alphaL)
    int_array[i, , , 4, 3, N] <- perc_boot_interval_grid(Mixture, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 5, 3, N] <- perc_boot_interval_grid(Student, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 6, 3, N] <- perc_boot_interval_grid(Triangular, parm = 2/3, B = bsim, alpha_length = alphaL)
    int_array[i, , , 7, 3, N] <- perc_boot_interval_grid(Uniform, parm = 1/2, B = bsim, alpha_length = alphaL)
    
    # Bias-corrected accelerated bootstrap
    int_array[i, , , 1, 4, N] <- bca_boot_interval_grid(Exponential, parm = 1, B = bsim, alpha_length = alphaL)
    int_array[i, , , 2, 4, N] <- bca_boot_interval_grid(Gaussian, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 3, 4, N] <- bca_boot_interval_grid(Lognormal, parm = 1.63, B = bsim, alpha_length = alphaL)
    int_array[i, , , 4, 4, N] <- bca_boot_interval_grid(Mixture, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 5, 4, N] <- bca_boot_interval_grid(Student, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 6, 4, N] <- bca_boot_interval_grid(Triangular, parm = 2/3, B = bsim, alpha_length = alphaL)
    int_array[i, , , 7, 4, N] <- bca_boot_interval_grid(Uniform, parm = 1/2, B = bsim, alpha_length = alphaL)
    
    # Bootstrap t
    int_array[i, , , 1, 5, N] <- bootstrap_t_grid(Exponential, parm = 1, B = bsim, alpha_length = alphaL)
    int_array[i, , , 2, 5, N] <- bootstrap_t_grid(Gaussian, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 3, 5, N] <- bootstrap_t_grid(Lognormal, parm = 1.63, B = bsim, alpha_length = alphaL)
    int_array[i, , , 4, 5, N] <- bootstrap_t_grid(Mixture, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 5, 5, N] <- bootstrap_t_grid(Student, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 6, 5, N] <- bootstrap_t_grid(Triangular, parm = 2/3, B = bsim, alpha_length = alphaL)
    int_array[i, , , 7, 5, N] <- bootstrap_t_grid(Uniform, parm = 1/2, B = bsim, alpha_length = alphaL)
    
    # Short bootstrap t
    int_array[i, , , 1, 6, N] <- short_bootstrap_t_grid(Exponential, parm = 1, B = bsim, alpha_length = alphaL)
    int_array[i, , , 2, 6, N] <- short_bootstrap_t_grid(Gaussian, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 3, 6, N] <- short_bootstrap_t_grid(Lognormal, parm = 1.63, B = bsim, alpha_length = alphaL)
    int_array[i, , , 4, 6, N] <- short_bootstrap_t_grid(Mixture, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 5, 6, N] <- short_bootstrap_t_grid(Student, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 6, 6, N] <- short_bootstrap_t_grid(Triangular, parm = 2/3, B = bsim, alpha_length = alphaL)
    int_array[i, , , 7, 6, N] <- short_bootstrap_t_grid(Uniform, parm = 1/2, B = bsim, alpha_length = alphaL)
    
    # Short percentile bootstrap
    int_array[i, , , 1, 7, N] <- short_percentile_bootstrap_grid(Exponential, parm = 1, B = bsim, alpha_length = alphaL)
    int_array[i, , , 2, 7, N] <- short_percentile_bootstrap_grid(Gaussian, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 3, 7, N] <- short_percentile_bootstrap_grid(Lognormal, parm = 1.63, B = bsim, alpha_length = alphaL)
    int_array[i, , , 4, 7, N] <- short_percentile_bootstrap_grid(Mixture, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 5, 7, N] <- short_percentile_bootstrap_grid(Student, parm = 0, B = bsim, alpha_length = alphaL)
    int_array[i, , , 6, 7, N] <- short_percentile_bootstrap_grid(Triangular, parm = 2/3, B = bsim, alpha_length = alphaL)
    int_array[i, , , 7, 7, N] <- short_percentile_bootstrap_grid(Uniform, parm = 1/2, B = bsim, alpha_length = alphaL)
    
    # Empirical Likelihood (Chi-square)
    int_array[i, , , 1, 8, N] <- emp_lik_grid(Exponential, parm = 1, type = "C", alpha_length = alphaL)
    int_array[i, , , 2, 8, N] <- emp_lik_grid(Gaussian, parm = 0, type = "C", alpha_length = alphaL)
    int_array[i, , , 3, 8, N] <- emp_lik_grid(Lognormal, parm = 1.63, type = "C", alpha_length = alphaL)
    int_array[i, , , 4, 8, N] <- emp_lik_grid(Mixture, parm = 0, type = "C", alpha_length = alphaL)
    int_array[i, , , 5, 8, N] <- emp_lik_grid(Student, parm = 0, type = "C", alpha_length = alphaL)
    int_array[i, , , 6, 8, N] <- emp_lik_grid(Triangular, parm = 2/3, type = "C", alpha_length = alphaL)
    int_array[i, , , 7, 8, N] <- emp_lik_grid(Uniform, parm = 1/2, type = "C", alpha_length = alphaL)
    
    # Empirical Likelihood (F)
    int_array[i, , , 1, 9, N] <- emp_lik_grid(Exponential, parm = 1, type = "F", alpha_length = alphaL)
    int_array[i, , , 2, 9, N] <- emp_lik_grid(Gaussian, parm = 0, type = "F", alpha_length = alphaL)
    int_array[i, , , 3, 9, N] <- emp_lik_grid(Lognormal, parm = 1.63, type = "F", alpha_length = alphaL)
    int_array[i, , , 4, 9, N] <- emp_lik_grid(Mixture, parm = 0, type = "F", alpha_length = alphaL)
    int_array[i, , , 5, 9, N] <- emp_lik_grid(Student, parm = 0, type = "F", alpha_length = alphaL)
    int_array[i, , , 6, 9, N] <- emp_lik_grid(Triangular, parm = 2/3, type = "F", alpha_length = alphaL)
    int_array[i, , , 7, 9, N] <- emp_lik_grid(Uniform, parm = 1/2, type = "F", alpha_length = alphaL)
    
    # Progress
    if (i %% (nsim / 5) == 0) {
      cat(paste0(i * 100 / nsim, "%"), "done with sample size", n, "simulations", "\n")
    }
  }
  # Increasing sample size indexing counter
  N <- N + 1
}
# Time check
end <- proc.time()
sim_time <- end - start
if(sim_time[3] < 60) {
  cat("Simulations took", sim_time[3], "seconds... \n")
} else if (sim_time[3] >= 60 && sim_time[3] < 3600) {
  cat("Simulations took", sim_time[3]/60, "minutes... \n")
} else {
  cat("Simulations took", sim_time[3]/3600, "hours... \n")
}

# Naming
dimnames(int_array) <- list(paste("Sim", 1:nsim),
                            paste("Alpha", seq(0.0005, 0.20, length = alphaL)),
                            c("Lower", "Upper", "Cap?", "Length", "Alpha"),
                            c("Exponential", "Gaussian", "Lognormal",
                              "Mixture", "Student's t", "Triangular", "Uniform"),
                            c("z-interval", "t-interval", "Percentile",
                              "BCa", "t", "Short t", "Short Percentile",
                              "Empirical (Chi)", "Empirical (F)"),
                            paste0("SS", nvec))

# Save object
# save(int_array, file = "./Resampling Code//SSS1.RData")
# save(int_array, file = "./Resampling Code//SSS2.RData")
# save(int_array, file = "./Resampling Code//SSS3.RData")
