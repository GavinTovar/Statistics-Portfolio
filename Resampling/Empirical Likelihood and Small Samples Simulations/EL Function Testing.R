### Empirical Likelihood Confidence Interval for a Univariate Mean
library(Rcpp)
sourceCpp("./Resampling Code/Small_Sims_Code.cpp")
# Weight function
w_i <- function(X, gamma, i) {
  a_i <- 1 / (X[i] + gamma)
  S <- sum(1 / (X + gamma))
  return(a_i / S)
}

# Vectorized version of above function
w_vec <- function(X, gamma){
  return(((X + gamma)*sum(1 / (X + gamma)))^(-1))
}

# First derivative of weight
w_i_prime <- function(X, gamma, i) {
  a_i <- 1 / (X[i] + gamma)
  a_i_prime <- 1 / (X[i] + gamma)^2
  S <- sum(1 / (X + gamma))
  S_prime <- sum(1 / (X + gamma)^2)
  - a_i_prime / S + a_i * S_prime / S^2
}

# Vectorized version of above function
w_vec_prime <- function(X, gamma){
  a_term <- 1 / (X + gamma)
  a_term_prime <- 1 / (X + gamma)^2
  S <- sum(1 / (X + gamma))
  S_term_prime <- sum(1 / (X + gamma)^2) 
  return(-a_term_prime / S + a_term * S_term_prime / S^2)
}

# Check domain validity 
valid_gamma <- function(X, gamma) {
  all(X + gamma > 0) || all(X + gamma < 0)
}

# f(gamma)
f_fun <- function(X, gamma, r0) {
  # Gamma check
  if (!valid_gamma(X, gamma)) return(NA_real_)
  
  n <- length(X)
  num <- 0
  for(i in 1:n) {
    num <- num + log(n * w_i(X, gamma, i))
  }
  return(-2 * num - r0)
}

# Vectorized version of above function
f_fun_vec <- function(X, gamma, r0){
  if (!valid_gamma(X, gamma)) return(NA_real_)
  
  n <- length(X)
  return(-2 * sum(log(n * w_vec(X, gamma))) - r0)
}


# f'(gamma)
f_prime_fun <- function(X, gamma, r0) {
  # Gamma check
  if (!valid_gamma(X, gamma)) return(NA_real_)
  
  n <- length(X)
  num <- 0
  for(i in 1:n) {
    num <- num + 1 / w_i(X, gamma, i) * w_i_prime(X, gamma, i)
  }
  return(-2 * num)
}

f_prime_fun_vec <- function(X, gamma, r0) {
  # Gamma check
  if (!valid_gamma(X, gamma)) return(NA_real_)
  
  return(-2 * sum(w_vec_prime(X, gamma) / w_vec(X, gamma)))
}

# Newton solver
Newton_Fun <- function(f, f_prime, gamma_start, X, r0, epsilon = 1e-7, max.iter = 100){
  gamma_new.vec <- rep(NA_real_, length(gamma_start))
  
  for (i in seq_along(gamma_start)) {
    iter <- 0
    gamma_old <- gamma_start[i]
    
    # Check for valid theta start
    if (gamma_start[i] < -max(X)) {
      gamma_new <- gamma_old - 100 * epsilon
      
    } else if (gamma_start[i] > -min(X)) {
      gamma_new <- gamma_old + 100 * epsilon
      
    } else {
      cat("Invalid choice of gamma_start; must be < -max(X) or > -min(X)")
      gamma_new <- NA_real_
      break
    }
    
    # print(paste("gamma_old = ", round(gamma_old, 6), ", gamma_new = ", round(gamma_new, 6), sep=""))
    while (!is.na(gamma_new) && !is.infinite(gamma_new) && abs(gamma_new - gamma_old) >= epsilon) {
      
      gamma_old <- gamma_new
      
      fval <- f(X, gamma_old, r0)
      fprime <- f_prime(X, gamma_old, r0)
      
      # print(paste("fval = ", round(fval, 3), ", fprime = ", round(fprime, 3), sep=""))
      
      
      if (is.na(fval) || is.na(fprime)) {
        cat("Invalid function value at iteration", iter, "for start", i, "\n")
        gamma_new <- NA_real_
        break
      }
      
      gamma_new <- gamma_old - fval / fprime
      # print(paste("gamma_old = ", round(gamma_old, 6), ", gamma_new = ", round(gamma_new, 6), sep=""))
      
      iter <- iter + 1
      if (iter > max.iter) {
        cat("Max iterations reached for starting value", i, "\n")
        gamma_new <- NA_real_
        break
      }
    }
    
    gamma_new.vec[i] <- gamma_new
  }
  return(data.frame(Starting_Val = gamma_start, Estimated_Theta = gamma_new.vec))
}

#####
# Simulations to explore CI performance
####

## Data-generating parameters
n <- 100
d <- 1

## Hypothesis Testing
nsim <- 1000
alpha <- 0.05
mu0_vec <- 1:10
r0 <- qchisq(1 - alpha, df = d)

## True parameter value
mu0 <- 4


####
# R function testing
####
muCIs <- matrix(0, nrow = nsim, ncol = 2)
for (i in 1:nsim) {
  # Sample
  samp <- rgamma(n, shape = 2, scale = 2)
  
  # Starting values for gamma
  starts <- c(-max(samp) - 0.01 * diff(range(samp)),
              -min(samp) + 0.01 * diff(range(samp)))
  
  gammaVals <- Newton_Fun(f_fun_vec, f_prime_fun_vec, starts, X = samp, r0 = r0, max.iter = 100)
  
  # Once gammas are obtained, plug back in to identify mu+ and mu-
  muCIs[i, 2] <- sum(w_vec(samp, gammaVals$Estimated_Theta[1]) * samp)
  muCIs[i, 1] <- sum(w_vec(samp, gammaVals$Estimated_Theta[2]) * samp)
}

# Coverage Rate
mean(muCIs[,1] < mu0 & mu0 < muCIs[,2])


####
# Cpp function testing
####
muCIs <- matrix(0, nrow = nsim, ncol = 2)
for (i in 1:nsim) {
  # Sample
  samp <- rgamma(n, shape = 2, scale = 2)
  
  # SCE: Need to modify these to make sure they are valid *gamma*
  # values. Need gamma to be less than -X_{(n)}  or greater than -X_{(1)}
  starts <- c(-max(samp) - 0.01 * diff(range(samp)),
              -min(samp) + 0.01 * diff(range(samp)))
  
  gammaVals <- Newton_Fun_cpp(f_fun_cpp, f_prime_fun_cpp, starts, X = samp, r0 = r0)
  
  # Once gammas are obtained, plug back in to identify mu+ and mu-
  muCIs[i, 2] <- sum(w_cpp(samp, gammaVals$Estimated_gamma[1]) * samp)
  muCIs[i, 1] <- sum(w_cpp(samp, gammaVals$Estimated_gamma[2]) * samp)
}

# Coverage Rate
mean(muCIs[, 1] < mu0 & mu0 < muCIs[, 2])

# Checking empirical likelihood function works
# emp_lik(samp, mu0, type = "C")
# emp_lik(samp, mu0, type = "F")
# emp_lik_grid(samp, mu0, type = "F")


####
# Testing cpp and R function equivalencies
####
samp <- rgamma(n, shape = 2, scale = 2)
starts <- c(-max(samp) - 0.01 * diff(range(samp)), -min(samp) + 0.01 * diff(range(samp)))
r0 <- qchisq(1 - alpha, df = d)

# w_cpp / w_vec
all.equal(w_vec(samp, starts[1]),
          w_cpp(samp, starts[1]))

# w_vec_prime / w_prime_cpp
all.equal(w_vec_prime(samp, starts[1]),
          w_prime_cpp(samp, starts[1]))

# f_fun_vec / f_fun_cpp
all.equal(f_fun_vec(samp, starts[1], r0),
          f_fun_cpp(samp, starts[1], r0))

# f_prime_fun_vec / f_prime_fun_cpp
all.equal(f_prime_fun_vec(samp, starts[1], r0),
          f_prime_fun_cpp(samp, starts[1], r0))

# Newton_Fun / Newton_Fun_cpp
all.equal(Newton_Fun(f_fun_vec, f_prime_fun_vec, starts, X = samp, r0 = r0, max.iter = 100),
          Newton_Fun_cpp(f_fun_cpp, f_prime_fun_cpp, starts, X = samp, r0 = r0, max_iter = 100))

#####
# Compare vectorized (between Cpp and R) to unvectorized versions in R:
#####

## Data-generating parameters
n <- 100
d <- 1

## Hypothesis Testing
nsim <- 1000
alpha <- 0.05
r0 <- qchisq(1 - alpha, df = d)

## Generate data
samp <- rgamma(n, shape = 2, scale = 2)

## Identify safe starts
starts <- c(-max(samp) - 0.01 * diff(range(samp)), -min(samp) + 0.01 * diff(range(samp)))

all.equal(f_fun(samp, starts[1], r0), f_fun_vec(samp, starts[1], r0))
all.equal(f_prime_fun(samp, starts[1], r0), f_prime_fun_vec(samp, starts[1], r0))

library(microbenchmark)
microbenchmark(f_fun(samp, starts[1], r0), 
               f_fun_vec(samp, starts[1], r0),
               f_fun_cpp(samp, starts[1], r0),
               times = 10000)
microbenchmark(f_prime_fun(samp, starts[1], r0),
               f_prime_fun_vec(samp, starts[1], r0),
               f_prime_fun_cpp(samp, starts[1], r0),
               times = 10000)
microbenchmark(Newton_R = Newton_Fun(f_fun_vec, f_prime_fun_vec, starts, X = samp, r0 = r0, max.iter = 100),
               Newton_Cpp = Newton_Fun_cpp(f_fun_cpp, f_prime_fun_cpp, starts, X = samp, r0 = r0, max_iter = 100),
               times = 1000)

####
# Testing cpp and R function equivalencies
####
samp <- rgamma(n, shape = 2, scale = 2)
starts <- c(-max(samp) - 0.01 * diff(range(samp)), -min(samp) + 0.01 * diff(range(samp)))
r0 <- qchisq(1 - alpha, df = d)

# w_cpp / w_vec
all.equal(w_vec(samp, starts[1]),
          w_cpp(samp, starts[1]))

# w_vec_prime / w_prime_cpp
all.equal(w_vec_prime(samp, starts[1]),
          w_prime_cpp(samp, starts[1]))

# f_fun_vec / f_fun_cpp
all.equal(f_fun_vec(samp, starts[1], r0),
          f_fun_cpp(samp, starts[1], r0))

# f_prime_fun_vec / f_prime_fun_cpp
all.equal(f_prime_fun_vec(samp, starts[1], r0),
          f_prime_fun_cpp(samp, starts[1], r0))

Newton_Fun_cpp(f_fun_cpp, f_prime_fun_cpp, starts, samp, r0) 
Newton_Fun(f_fun_vec, f_prime_fun_vec, starts, X = samp, r0 = r0, max.iter = 10)
