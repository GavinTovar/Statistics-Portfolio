### Set seed for reproducibility and load required packages
set.seed(10262025)
library(tidyverse)
library(mvtnorm)
library(lme4)
library(patchwork)

### Data Generation Function
Mixed.Model.DataGen <- function(
    FixEff = 4,                         # Number of fixed effects (including intercept)
    RanEff = 2,                         # Number of random effects
    Levels = rep(8, RanEff),            # Number of levels for each random effect (vector)
    N = 120,                            # Total sample size 
    Sigma.int = rep(3, RanEff),         # Std dev for random intercepts
    Sigma.slope = rep(1, RanEff),       # Std dev for random slopes 
    Cov.mat = NULL,                     # Covariance matrix for fixed-effect covariates
    Beta = NULL,                        # Fixed effect coefficients (vector)
    mu.vec = NULL,                      # Mean vector for covariates
    Response.distribution = "Normal",   # Error distribution type: "Normal", "Poisson", "Bernoulli", etc.
    Error.params = list(sd = 1),        # Parameters for error distribution (e.g., sd for normal, updated to 1)
    Link.function = "identity",         # Link function for GLMM (e.g., "identity", "log", "logit")
    overdispersion = FALSE              # Attempt to induce overdispersion through hierarchical methods
) {
  # Input checks
  RanEff <- length(Levels)
  if (length(Sigma.int) != RanEff || length(Sigma.slope) != RanEff) {
    stop("Length of Sigma.int and Sigma.slope must match the number of random effects (RanEff).")
  }
  
  # Calculate sample size per level by dividing total sample size N by the number of levels in each factor
  Level.SS <- sapply(Levels, function(levels) N / levels)
  if (any(Level.SS != round(Level.SS))) {
    stop("The total sample size (N) cannot be equally divided by the number of levels. Adjust N or Levels.")
  }
  
  # Default settings
  num.cov <- FixEff - 1  # Exclude intercept
  
  # Mean vector
  if (is.null(mu.vec)) {
    mu.vec <- rep(0, num.cov)
  }
  if (length(mu.vec) != num.cov) {
    stop("mu.vec must have length equal to the number of fixed-effect covariates (FixEff - 1).")
  }
  
  # Covariance matrix
  if (is.null(Cov.mat)) {
    Cov.mat <- matrix(0.5, nrow = length(mu.vec), ncol = length(mu.vec))
    diag(Cov.mat) <- 1
  }
  if (!all(dim(Cov.mat) == num.cov)) {
    stop("Cov.mat dimensions must match the number of fixed-effect covariates (FixEff - 1).")
  }
  
  # Fixed effect beta sequence
  if (is.null(Beta)) {
    Beta <- c(1, rep(1, num.cov))
  }
  if (length(Beta) != FixEff) {
    stop("Length of Beta must equal the number of fixed effects (FixEff).")
  }
  
  # Generate fixed-effect covariates
  x.mat <- matrix(rnorm(N * num.cov, mean = mu.vec), nrow = N, ncol = num.cov)
  colnames(x.mat) <- paste0("X", 1:num.cov)
  
  # Generate random effect levels and factors
  Factors <- lapply(seq_len(RanEff), function(i) {
    rep(factor(rep(seq_len(Levels[i]), each = Level.SS[i])), length.out = N)
  })
  names(Factors) <- paste0("Z", 1:RanEff)
  
  # Generate random intercepts and slopes
  Random.Eff <- lapply(seq_len(RanEff), function(i) {
    intercepts <- rnorm(Levels[i], mean = 0, sd = Sigma.int[i])
    slopes <- if (Sigma.slope[i] > 0) {
      rnorm(Levels[i], mean = 0, sd = Sigma.slope[i])
    } else {
      rep(0, Levels[i])
    }
    list(intercepts = intercepts[as.integer(Factors[[i]])], slopes = slopes[as.integer(Factors[[i]])])
  })
  
  
  #####
  
  ### Perhaps unneeded randomization element???
  # Randomize order of random effect variables and align Factors
  Random.Eff <- lapply(seq_along(Random.Eff), function(i) {
    # Randomize order of intercepts and slopes
    random_indices <- sample(seq_along(Random.Eff[[i]]$intercepts))
    list(
      intercepts = Random.Eff[[i]]$intercepts[random_indices],
      slopes = Random.Eff[[i]]$slopes[random_indices],
      random_indices = random_indices # Store the order for aligning Factors
    )
  })
  
  # Update Factors to match the randomized order
  Factors <- lapply(seq_along(Factors), function(i) {
    factor(as.integer(Factors[[i]])[Random.Eff[[i]]$random_indices])
  })
  names(Factors) <- paste0("Z", seq_along(Factors))
  
  #####
  
  
  # Construct the linear predictor (fixed effects + random effects)
  lin.pred <- Beta[1] + x.mat %*% Beta[-1]
  for (i in seq_len(RanEff)) {
    lin.pred <- lin.pred + Random.Eff[[i]]$intercepts + Random.Eff[[i]]$slopes * x.mat[, 1]
  }
  
  # Apply link function (if GLMM)
  if (Link.function == "log") {
    lin.pred <- exp(lin.pred)
    if(any(lin.pred > 100)){
      print("Woah there bucko -- Inputs are creating very large counts, try making them smaller if this is not intended")
    }
  } else if (Link.function == "logit") {
    lin.pred <- 1 / (1 + exp(lin.pred))
    if(any(lin.pred > 0.95)){
      print("Woah there bucko -- Inputs are creating very large success probabilities, try making them smaller if this is not intended")
    }
  }
  
  # Construct the response variable
  y.vec <- switch(Response.distribution,
                  Normal = rnorm(N, mean = lin.pred, sd = Error.params$sd),
                  Poisson = rpois(N, lambda = lin.pred),
                  Bernoulli = rbinom(N, size = 1, prob = lin.pred),
                  stop("Unsupported error distribution")
  )
  
  # Poisson overdispersion attempt
  # Thinking of the current Y values as groups and then each group has counts distributed poisson
  if(Response.distribution == "Poisson" & overdispersion == TRUE){
    z.vec <- rep(0, N)
    for(i in 1:N){
      z.vec[i] <- sum(rpois(y.vec[i], 3)) # Lambda value here is completely arbitrary
    }
    y.vec <- z.vec
  }
  
  # Combine data into a dataframe
  df <- data.frame(Y = y.vec, x.mat, do.call(cbind, Factors))
  return(df)
}

### Model Selection Testing
library(MuMIn)
nsim <- 1000

# Model selection function
Post.Model.GLMM_Coeffs <- function(full.model) {
  
  # Perform model selection using dredge
  model_selection <- dredge(full.model, beta = "none")
  
  # Get models ranked by AICc
  mod.list <- get.models(model_selection, subset = TRUE)
  
  # Pull out the rank 1 model
  Select.Model <- summary(eval(summary(mod.list[[1]])$call))
  
  # Extract fixed effect coefficients
  Fixed.Coeff <- Select.Model$coefficients
  
  # Extract random effect coefficients
  Random.Coeff <- unlist(Select.Model$varcor)
  
  # Return results
  return(list(Fixed.Coeff = Fixed.Coeff, Random.Coeff = Random.Coeff))
}

# for(i in 1:nsim)
sim.data <- Mixed.Model.DataGen(
  FixEff = 2,
  RanEff = 1,
  Levels = c(3),
  N = 120,
  Sigma.int = c(0.5),
  Sigma.slope = c(0),
  Cov.mat = NULL,
  Beta = c(0.25, 0.25),
  mu.vec = c(1),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE
)
# Fit GLMM
PoisMod <- glmer(Y ~ X1 + (X1 | Z1), 
                 data = sim.data, 
                 family = poisson(link = "log"),
                 na.action = na.fail)

# Perform model selection 
Model_Coeffs <- Post.Model.GLMM_Coeffs(PoisMod)

# Correct Model
CorrectPoisMod <- glmer(Y ~ X1 + (1 | Z1), 
                 data = sim.data, 
                 family = poisson(link = "log"))

summary(CorrectPoisMod)