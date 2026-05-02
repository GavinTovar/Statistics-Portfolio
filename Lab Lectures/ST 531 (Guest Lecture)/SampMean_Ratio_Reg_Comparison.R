# Comparison of Sample mean, Ratio Estimator, and Regression Estimator

# Set seed for reproducibility
set.seed(12345)

# Population size
N <- 1000

# Sample size
n <- 30

#####
# Generate Population: Run once
#####

# Generate auxiliary variable x
mean_x <- 10
sd_x <- 1
x_pop <- rnorm(N, mean_x, sd_x)
x_pop_mean <- mean(x_pop)

# Generate random errors
sd_sig <- 1
errs <- rnorm(N, 0, sd_sig)

# Generate target variable y
beta_0 <- 5
beta_1 <- 0.5
y_pop <- beta_0 + beta_1 * x_pop + errs
y_pop_mean <- mean(y_pop)

#####
# Explore sampling behavior of estimators
#####

# Set number of samples to draw to approximate the sampling distribution
nsamp <- 100000

# Create empty matrices to store the estimates (sample mean, ratio,
# regression) and the standard errors of those estimates
ests <- matrix(0, nrow = nsamp, ncol = 3)
se_ests <- matrix(0, nrow = nsamp, ncol = 3)
colnames(ests) <- c("SampleMean", "Ratio", "Regression")
colnames(se_ests) <- c("SampleMean", "Ratio", "Regression")

# Iterate through drawing new samples,
# calculating estimators and SEs for each
for (i in 1:nsamp) {
  sampInd <- sample(1:N, n, replace = FALSE)
  
  x_samp <- x_pop[sampInd]
  y_samp <- y_pop[sampInd]
  Bhat <- mean(y_samp) / mean(x_samp)
  Bhat1 <- cov(x_samp, y_samp) / var(x_samp)
  Bhat0 <- mean(y_samp) - Bhat1 * mean(x_samp)
  
  ests[i, 1] <- mean(y_samp)
  ests[i, 2] <- Bhat * x_pop_mean
  ests[i, 3] <- Bhat0 + Bhat1 * x_pop_mean
  
  resids_rat <- y_samp - Bhat * x_samp
  resids_reg <- y_samp - (Bhat0 + Bhat1 * x_samp)
  
  s2_resids_rat <- var(resids_rat)
  s2_resids_reg <- 1 / (n - 2) * sum((resids_reg - mean(resids_reg))^2)
  
  fpc <- (1 - n / N)
  se_ests[i, 1] <- sqrt(fpc * var(y_samp) / n)
  se_ests[i, 2] <- sqrt(fpc * (x_pop_mean / mean(x_samp))^2 * s2_resids_rat / n)
  se_ests[i, 3] <- sqrt(fpc * s2_resids_reg / n)
}

# Mean of estimators
colMeans(ests)

# Bias of estimators
colMeans(ests - y_pop_mean)

# Mean squared error of estimators
colMeans((ests - y_pop_mean)^2)

# Mean of standard errors
colMeans(se_ests)

# Sampling standard deviation of estimators
apply(ests, 2, sd)

# Confidence interval coverage
alpha <- 0.05
crit.val <- qnorm(1 - alpha / 2)

CI_LB <- ests - crit.val * se_ests
CI_UB <- ests + crit.val * se_ests

apply((CI_LB < y_pop_mean & CI_UB > y_pop_mean), 2, mean)