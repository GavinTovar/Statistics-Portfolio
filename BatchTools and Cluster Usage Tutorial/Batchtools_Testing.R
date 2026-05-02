### Batchtools Testing Script
library(batchtools)
removeRegistry(wait = 5)
reg <- makeExperimentRegistry(file.dir = "./BatchtoolsTestingRegistry", conf.file = ".batchtools.conf.R")
# reg <- loadRegistry("./BatchtoolsTestingRegistry", writeable = TRUE)

# ---------------------------
# Problem: Data generating process
# ---------------------------
addProblem(
  name = "MLRdata",
  data = list(Beta = c(1, 0.5, -0.25, 0.75)),
  fun = function(job, data, n, rho, error_variance) {
    # Ensure mvtnorm is available
    if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Need mvtnorm package")
    
    # Instantiate variables that are constant over all jobs
    Beta <- data$Beta
    
    # Create correlation matrix
    p = length(Beta) - 1
    Sigma <- matrix(rho, nrow = p, ncol = p)
    diag(Sigma) <- 1
    
    # Generate predictors ~ MVN(0, Sigma)
    X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    colnames(X) <- paste0("X", 1:p)
    
    # Create response vector
    y <- Beta[1] + X %*% Beta[-1] + rnorm(n, sd = sqrt(error_variance))
    
    return(data.frame(y = y, X))
  }
)

# ---------------------------
# Algorithm 1: OLS
# ---------------------------
addAlgorithm(
  name = "ols",
  fun = function(job, data, instance) {
    fit <- lm(y ~ ., data = instance)
    coefs <- coef(summary(fit))
    list(
      method   = "ols",
      beta_hat = coefs[, "Estimate"],
      se_hat   = coefs[, "Std. Error"],
      rse      = summary(fit)$sigma
    )
  }
)

# ---------------------------
# Algorithm 2: Ridge regression
# ---------------------------
addAlgorithm(
  name = "ridge",
  fun = function(job, data, instance, lambda) {
    # Ensure glmnet is available
    if (!requireNamespace("glmnet", quietly = TRUE)) stop("Need glmnet package")
    
    X <- as.matrix(instance[, -1]) # predictors
    y <- instance$y
    
    # pars$lambda comes from algo.designs
    fit <- glmnet::glmnet(X, y, alpha = 0, lambda = lambda)
    
    list(
      method   = "ridge",
      lambda   = lambda,
      beta_hat = as.numeric(coef(fit))
    )
  }
)

# ---------------------------
# Experimental designs
# ---------------------------
prob.designs <- list(
  MLRdata = expand.grid(
    n = c(50, 100),
    rho = c(0.3, 0.6),
    error_variance = c(1, 4)
  )
)

algo.designs <- list(
  ols = data.frame(),  # no tuning parameters
  ridge = expand.grid(lambda = c(0.01, 0.1, 1))  # tuning parameter grid
)


# Add jobs = all combinations
addExperiments(
  prob.designs,
  algo.designs
)

# Run jobs
submitJobs()

# While waiting on jobs:
# getStatus()

# If any error
# getErrorMessages()

# Breakdown of problem-experiment combinations
# summarizeExperiments()

# When jobs have finished:
results <- reduceResultsList()



####### Batchmap() approach:
reg <- makeRegistry(file.dir = "./BatchtoolsTestingRegistry", conf.file = ".batchtools.conf.R")

# ---------------------------
# Define function for one job
# ---------------------------
sim_fun <- function(n, rho, error_variance, method, lambda = NA) {
  # Ensure needed packages are installed
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Need mvtnorm package")
  }
  if (method == "ridge" && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Need glmnet package")
  }
  # True coefficients
  Beta <- c(1, 0.5, -0.25, 0.75)
  p <- length(Beta) - 1
  
  # Correlation matrix
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  # Predictors ~ MVN(0, Sigma)
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  
  # Response
  y <- Beta[1] + X %*% Beta[-1] + rnorm(n, sd = sqrt(error_variance))
  dat <- data.frame(y = y, X)
  
  # Fit model depending on method
  if (method == "ols") {
    fit <- lm(y ~ ., data = dat)
    coefs <- coef(summary(fit))
    return(list(
      method   = "ols",
      beta_hat = coefs[, "Estimate"],
      se_hat   = coefs[, "Std. Error"],
      rse      = summary(fit)$sigma
    ))
  } else if (method == "ridge") {
    fit <- glmnet::glmnet(as.matrix(dat[, -1]), dat$y, alpha = 0, lambda = lambda)
    return(list(
      method   = "ridge",
      lambda   = lambda,
      beta_hat = as.numeric(coef(fit))
    ))
  }
}

# ---------------------------
# Parameter grid
# ---------------------------
param_grid <- expand.grid(
  n = c(50, 100),
  rho = c(0.3, 0.6),
  error_variance = c(1, 4),
  method = c("ols", "ridge"),
  lambda = c(NA, 0.01, 0.1, 1),
  stringsAsFactors = FALSE
)

# Filter out invalid combos: only ridge uses lambda
param_grid <- subset(param_grid, (method == "ols" & is.na(lambda)) |
                       (method == "ridge" & !is.na(lambda)))

# ---------------------------
# Map parameters to jobs
# ---------------------------
batchMap(fun = sim_fun, 
         n = param_grid$n, 
         rho = param_grid$rho,
         error_variance = param_grid$error_variance,
         method = param_grid$method,
         lambda = param_grid$lambda)

# Submit jobs
submitJobs()
