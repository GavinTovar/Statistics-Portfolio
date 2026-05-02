### Empirical Likelihood Inference for a Univariate Mean
## Data-generating parameters
n <- 100
d <- 1

## Hypothesis Testing
nsim <- 1000
alpha <- 0.05
mu0_vec <- 1:10
r0 <- qchisq(1 - alpha, df = d)
error_count <- 0
reject_mat <- matrix(0, nrow = nsim, ncol = length(mu0_vec))

# Solution for weights to maximize log of empirical likelihood ratio
wtsFun <- function(samp, lam, mu0) {
  n <- length(samp)
  wts <- (1 / n) * (1 / (1 + lam * (samp - mu0)))
  return(wts)
}

# Number of lambda values to evaluate
numLam <- 1000
logelr_stat <- matrix(0, nrow = nsim, ncol = length(mu0_vec))
for (i in 1:nsim) {
  # Generate a new sample
  samp <- rgamma(n, shape = 2, scale = 2)
  for (j in 1:length(mu0_vec)) {
    # Compute the log empirical likelihood ratio
    if (mu0_vec[j] < min(samp) || mu0_vec[j] > max(samp)) {
      log_elr <- -Inf
      
    } else if (all(samp == mu0_vec[j])) {
      log_elr <- 0
      
    } else {
      # # Setting up grid given initial bounds
      lamMin <- (1 - (1 / n)) / (mu0_vec[j] - max(samp))
      lamMax <- (1 - (1 / n)) / (mu0_vec[j] - min(samp))
      lamGrid <- seq(lamMin, lamMax, length = numLam)
      
      # Fabricate estimating equation result vector
      estim_eq_vals <- numeric(numLam)
      
      # Calculate estimating equation values for each lambda
      for (k in 1:numLam) {
        denom <- 1 + lamGrid[k] * (samp - mu0_vec[j])
        
        # Check for non-positive denominators
        if(any(denom <= 0)) {
          estim_eq_vals[k] <- Inf
          next
        }
        # Calculate the estimating equation value
        estim_eq_vals[k] <- (1 / n) * sum((samp - mu0_vec[j]) / denom)
      }
      
      # Find the lambda that minimizes the absolute value of the estimating equation
      best_idx <- which.min(abs(estim_eq_vals))
      opt_lam <- lamGrid[best_idx]
      opt_wts <- wtsFun(samp, opt_lam, mu0_vec[j])
      
      # Check weight sum
      if (abs(sum(opt_wts) - 1) > 1e-6) {
        error_count <- error_count + 1
        warning("Weights do not sum to 1!")
      }
      # Compute the log empirical likelihood ratio
      log_elr <- sum(log(n * opt_wts))
      
    }
    # Add ELR scaling factor
    logelr_stat[i,j] <- -2 * log_elr
    
    # Reject null if log empirical likelihood ratio is greater than critical value
    reject_mat[i, j] <- logelr_stat[i, j] > r0
  }
  if (i %% (nsim / length(mu0_vec)) == 0) {
    cat(paste0(i * 100 / nsim, "% Done,"), "Simulation", i, "completed! \n")
  }
}

Results <- colMeans(reject_mat)
names(Results) <- paste0("mu0 = ", mu0_vec)
Results

## Confidence Intervals
WtsFun <- function(samp, gamma) {
  n <- length(samp)
  frac <- 1 / (samp + gamma)
  Wts <- frac / sum(frac)
  if (any(Wts < 0)) {
    return(NULL)
  } else {
    return(Wts)
  }
  
}

course_num <- 1000
fine_num <- course_num*5
mu <- 4
captured_course <- captured_fine <- logical(nsim)

for (i in 1:nsim) {
  # Generate a new sample
  samp <- rgamma(n, shape = 2, scale = 2)
  
  # Course grid search for mu- and mu+
  course_grid_minus <- seq(min(samp), by = -0.10, length = course_num)
  course_grid_plus <- seq(max(samp), by = 0.10, length = course_num)
  
  # Fabricate error objects
  mu_minus_error <- numeric(course_num)
  mu_plus_error <- numeric(course_num)
  for (j in 1:course_num) {
    mu_minus_error[j] <- abs(-2 * sum(log(n * WtsFun(samp, course_grid_minus[j]))) - log(r0))
    mu_plus_error[j] <- abs(-2 * sum(log(n * WtsFun(samp, course_grid_plus[j]))) - log(r0))
  }
  
  # Finding course min and max mu values
  start_mu_minus <- course_grid_minus[which.min(mu_minus_error)]
  start_mu_plus <- course_grid_plus[which.min(mu_plus_error)]
  
  # Check if mu is within the confidence interval
  captured_course[i] <- start_mu_minus < mu && mu < start_mu_plus
  
  # Fine grid search for mu- and mu+
  fine_grid_minus <- seq(start_mu_minus, by = -0.005, length = fine_num)
  fine_grid_plus <- seq(start_mu_plus, by = 0.005, length = fine_num)
  
  # Fabricate error objects
  mu_minus2_error <- numeric(fine_num)
  mu_plus2_error <- numeric(fine_num)
  for (k in 1:course_num) {
    mu_minus2_error[k] <- abs(-2 * sum(log(n * WtsFun(samp, fine_grid_minus[k]))) - log(r0))
    mu_plus2_error[k] <- abs(-2 * sum(log(n * WtsFun(samp, fine_grid_plus[k]))) - log(r0))
  }
  
  # Finding min and max mu values
  mu_minus <- fine_grid_minus[which.min(mu_minus2_error)]
  mu_plus <- fine_grid_plus[which.min(mu_plus2_error)]
  
  # Check if mu is within the confidence interval
  captured_fine[i] <- mu_minus < mu && mu < mu_plus
  
  # Print Progress
  if (i %% (nsim / 5) == 0) {
    cat(i * 100 / nsim, "% Done,", "Simulations", i, "completed! \n")
  }
}
mean(captured_course)
mean(captured_fine)

#########
iter_search <- function(samp, mu, r0, n_iter = 5, n_grid = 200, step_init = 0.2) {
  n <- length(samp)
  
  search_one_side <- function(start, direction) {
    step <- step_init
    best <- start
    for (iter in 1:n_iter) {
      grid <- seq(best, by = direction * step, length.out = n_grid)
      errs <- sapply(grid, function(g) {
        Wts <- WtsFun(samp, g)
        if (is.null(Wts)) return(Inf)
        abs(-2 * sum(log(n * Wts)) - log(r0))
      })
      best <- grid[which.min(errs)]
      step <- step / 5   # refine step
    }
    best
  }
  
  mu_minus <- search_one_side(min(samp), -1)
  mu_plus  <- search_one_side(max(samp), +1)
  
  list(mu_minus = mu_minus, mu_plus = mu_plus,
       captured = (mu_minus < mu && mu < mu_plus))
}
iter_search(samp, mu = 4, qchisq(1 - alpha, df = 1), n_iter = 10, n_grid = 100, step_init = 0.0001)

### The course grid is supposed to help find a starting pair for
### Newton's method to replace the fine grid search

# First derivative function of w_i w.r.t. gamma
w_i_prime <- function(X, gamma, i) {
  a_i <- 1 / (X[i] + gamma)
  a_i_prime <- -1 / (X[i] + gamma)^2
  S <- sum(1 / (X + gamma))
  S_prime <- sum(-1 / (X + gamma)^2)
  
  numerator <- a_i_prime * S - a_i * S_prime
  denominator <- S^2
  
  return(numerator / denominator)
}

# Second derivative function of w_i w.r.t. gamma
w_i_double_prime <- function(X, gamma, i) {
  xg <- X + gamma
  a_i <- 1 / xg[i]
  a_i_prime <- -1 / xg[i]^2
  a_i_double_prime <- 2 / xg[i]^3
  
  S <- sum(1 / xg)
  S_prime <- sum(-1 / xg^2)
  S_double_prime <- sum(2 / xg^3)
  
  numerator <- a_i_double_prime * S^2 - 2 * a_i_prime * S * S_prime -
    a_i * S_double_prime * S + 2 * a_i * S_prime^2
  denominator <- S^3
  
  return(numerator / denominator)
}


Newton.Fun <- function(d1, d2, theta.start, epsilon = 1e-7, max.iter = 100, obs.x){
  theta.new.vec <- rep(0, length(theta.start))
  for(i in 1:length(theta.start)){
    # print(paste("Starting Value", i, "iterations"))
    iter <- 0
    theta.old <- theta.start[i]
    theta.new <- theta.start[i] + 100 * epsilon
    
    # Newton-Raphson iterative updates
    while(abs(theta.new - theta.old) >= epsilon){
      theta.old <- theta.new
      theta.new <- theta.old - d1(theta.old, x = obs.x) / d2(theta.old, x = obs.x)  
      
      # Keeping track of iterations
      iter <- iter + 1
      if (iter > max.iter) {
        cat("Max iterations reached for starting value", i, "\n")
        # Assign Inf or -Inf to diverging values
        theta.new <- ifelse(theta.new > 0, Inf, -Inf)
        break
      }
      # print(theta.new)
    }
    theta.new.vec[i] <- theta.new
  }
  return(data.frame(Starting_Val = theta.start, Estimated_Theta = theta.new.vec))
}
ci_coverage <- logical(nsim)
# Function to compute log empirical likelihood ratio
logELR <- function(samp, mu0) {
  n <- length(samp)
  
  # Skip if out of range
  if (mu0 < min(samp) || mu0 > max(samp)) return(Inf)
  
  # Grid search over lambda
  lamMin <- (1 - 1/n) / (mu0 - max(samp))
  lamMax <- (1 - 1/n) / (mu0 - min(samp))
  lamGrid <- seq(lamMin, lamMax, length.out = 1000)
  
  est_eq_vals <- numeric(length(lamGrid))
  for (k in seq_along(lamGrid)) {
    denom <- 1 + lamGrid[k] * (samp - mu0)
    if (any(denom <= 0)) {
      est_eq_vals[k] <- Inf
    } else {
      est_eq_vals[k] <- mean((samp - mu0) / denom)
    }
  }
  
  best_idx <- which.min(abs(est_eq_vals))
  opt_lam <- lamGrid[best_idx]
  wts <- wtsFun(samp, opt_lam, mu0)
  
  if (any(wts <= 0)) return(Inf)
  
  # Log empirical likelihood ratio
  return(-2 * sum(log(n * wts)))
}

for (i in 1:nsim) {
  samp <- rgamma(n, shape = 2, scale = 2)
  
  # Search over a grid of mu values around sample mean
  mu_grid <- seq(2, 6, length.out = 500)
  elr_vals <- sapply(mu_grid, function(mu0) logELR(samp, mu0))
  
  # CI is where ELR <= critical value
  in_ci <- mu_grid[elr_vals <= r0]
  
  if (length(in_ci) > 0) {
    mu_lo <- min(in_ci)
    mu_hi <- max(in_ci)
    ci_coverage[i] <- (mu_lo <= mu) && (mu <= mu_hi)
  } else {
    ci_coverage[i] <- FALSE
  }
  
  if (i %% (nsim / 5) == 0) {
    cat(i * 100 / nsim, "% Done,", "Simulations", i, "completed!\n")
  }
}

mean(ci_coverage)