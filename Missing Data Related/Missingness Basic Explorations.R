### (ST 599) Homework 1
library(Surrogate)
data("Schizo_PANSS")
library(knitr)
library(kableExtra)
# Schizo_PANSS <- read.csv("Schizo_PANSS.csv")

hw_data <- Schizo_PANSS[,c("Id", "Treat", "Week1", "Week4", "Week8")]

### Part 1
# Identify unique missingness patterns
Miss.Mat <- unique(is.na(hw_data))

# Add a column for proportions, initialized to 0
Miss.Mat <- cbind.data.frame(Proportion = rep(0, nrow(Miss.Mat)), Miss.Mat)

# Calculate the proportion of each missingness pattern
for (i in seq_len(nrow(Miss.Mat))) {
  # Extract the pattern as a vector
  pattern <- as.logical(Miss.Mat[i, -1])
  
  # Compare rows of hw_data with the pattern
  matches <- apply(is.na(hw_data), 1, function(row) all(row == pattern))
  
  # Compute the proportion of rows matching this pattern
  Miss.Mat$Proportion[i] <- sum(matches) / nrow(hw_data)
}

Miss.Mat


# Create a formatted kable for the results
kable(Miss.Mat, format = "latex", booktabs = TRUE,
      col.names = c("Proportion", colnames(hw_data)),
      caption = "Missingness Patterns and Proportions in Schizo_PANSS")

### Part 2
# Creating a factor variable based on missingness pattern
Pattern.Match <- function(Vec){
  if(all(is.na(Vec) == Miss.Mat[1, -1])){
    return("FFFFF")
    
  } else if(all(is.na(Vec) == Miss.Mat[2, -1])){
    return("FFFFT")
    
  } else if(all(is.na(Vec) == Miss.Mat[3, -1])) {
    return("FFFTT")
    
  } else if(all(is.na(Vec) == Miss.Mat[4, -1])) {
    return("FFFTF")
    
  } else if(all(is.na(Vec) == Miss.Mat[5, -1])) {
    return("FFTFT")
    
  } else if(all(is.na(Vec) == Miss.Mat[6, -1])) {
    return("FFTTT")
    
  } else if(all(is.na(Vec) == Miss.Mat[7, -1])) {
    return("FFTFF")
    
  } else {
    return("FFTTF")
  }
}

HW.Data <- hw_data
HW.Data$Missingness <- apply(hw_data, 1, Pattern.Match)
HW.Data$Missingness <- factor(HW.Data$Missingness)
str(HW.Data)

# Multinominal Regression
library(nnet)
Multi.Mod <- multinom(Missingness ~ Treat, 
                      data = HW.Data, family = binomial)
summary(Multi.Mod)

Multi.Coef <- cbind.data.frame(Treat = summary(Multi.Mod)[["coefficients"]][, 1],
                               Std.Error = summary(Multi.Mod)[["standard.errors"]][, 2],
                                  z.stat = summary(Multi.Mod)[["coefficients"]][, 2] / 
                                    summary(Multi.Mod)[["standard.errors"]][, 2])
Multi.Coef <- cbind.data.frame(Multi.Coef, p.value = format(2*round(pnorm(abs(Multi.Coef[["z.stat"]]), lower.tail = FALSE), 6), 
                                                                  scientific = FALSE))
Multi.Coef

Multi.Coef %>%
  kable(
    format = "latex", 
    digits = 2, 
    align = "c", 
    booktabs = TRUE, 
    col.names = c("Treatment", "Std. Error", "Z-stat", "P-value"),
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"), # Only minimal styling
    position = "center"
  )

# Capture the summary of the model
mod_summary <- summary(Multi.Mod)

# Convert the summary coefficients into a data frame for better readability
coef_df <- as.data.frame(mod_summary$coefficients)

# Display the table using knitr
kable(coef_df, format = "latex", booktabs = TRUE,
      caption = "Multinomial Regression Coefficients for Missingness Patterns")

# Creating table of counts
Patterns <- unique(HW.Data$Missingnes)
Chi.Mat <- matrix(0, nrow = length(Patterns), ncol = 2)
for(i in seq_along(Patterns)){
  Chi.Mat[i, ] <-  c(nrow(subset(HW.Data, Treat == 1 & Missingness == Patterns[i])),
                     nrow(subset(HW.Data, Treat == -1 & Missingness == Patterns[i])))
}
colnames(Chi.Mat) <- c("Treatment", "No Treatment")
row.names(Chi.Mat) <- Patterns

# Get the expected counts from chi-squared test
expected_counts <- chisq.test(Chi.Mat)[["expected"]]

# Combine observed and expected counts into a single data frame
combined_table <- data.frame(
  Missingness = rownames(Chi.Mat),
  Treatment_Observed = Chi.Mat[, 1],
  No_Treatment_Observed = Chi.Mat[, 2],
  Treatment_Expected = round(expected_counts[, 1], 2),
  No_Treatment_Expected = round(expected_counts[, 2], 2)
)

# Create a kable output for side-by-side table
combined_table %>%
  kable(
    format = "latex", 
    col.names = c("Missingness Pattern", "Treatment (Observed)", "No Treatment (Observed)", "Treatment (Expected)", "No Treatment (Expected)"), 
    digits = 2, 
    align = "c", 
    booktabs = TRUE, 
    caption = "Observed and Expected Counts"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    position = "center"
  )

# Chi-Squared Test
chisq.test(Chi.Mat)

### Part 3
# How many people dropped out of the study after Week 1?
Miss.Mat[3, 1] * nrow(hw_data)

# How many people dropped out of the study after Week 4?
Miss.Mat[2, 1] * nrow(hw_data)

# For patients who dropped out at wk 4 is there evidence that 
# PANSS at the prior measurement predicted dropout?
Logistic.Data <- HW.Data
Logistic.Data[["Wk4Drop"]] <- ifelse(Logistic.Data[["Missingness"]] == Patterns[3], 1, 0)

# Logistic Regression to predict future missingness
summary(glm(Wk4Drop ~ Week1 + Treat, data = Logistic.Data))

# For patients who dropped out at wk 4 is there evidence that 
# PANSS at the prior measurement predicted dropout?
Logistic.Data[["Wk8Drop"]] <- ifelse(Logistic.Data[["Missingness"]] == Patterns[2], 1, 0)

# Logistic Regression to predict future missingness
summary(glm(Wk8Drop ~ Week4 + Treat, data = Logistic.Data))

# How many people had intermittent missingness?
sum(Miss.Mat[4:8, 1] * nrow(hw_data))

# For patients with intermittent missingness is there evidence that 
# PANSS at the prior measurement predicted dropout?
Logistic.Data[["IntDrop"]] <- ifelse(Logistic.Data[["Missingness"]] == Patterns[4] | Logistic.Data[["Missingness"]] == Patterns[5], 1, 0)

# Logistic Regression to predict future missingness
summary(glm(IntDrop ~ Week1 + Week8 + Treat, data = Logistic.Data))

### Part 4
# Read Rob's mind

### Part 5
## Algorithm construction and testing
# Test Data
set.seed(123)
n <- 10000
p <- 5
K <- 3 
X <- list()
beta <- rnorm(p)
L <- matrix(rnorm(K * K), K, K)
Sigma <- L %*% t(L)
y <- list()
for (i in 1:n) {
  X[[i]] <- matrix(rnorm(p * K), K, p)
  y[[i]] <- X[[i]] %*% beta + MASS::mvrnorm(mu = rep(0, K), Sigma = Sigma)
}

Algorithm <- function(X, y, Beta0, Sigma0, max.iter = 100, tol = 1e-11, verbose = FALSE){
  # Initialize parameters
  n <- length(y)
  Beta <- Beta0
  Sigma <- Sigma0
  
  # Initialize function objects
  iter <- 0
  converged1 <- FALSE
  converged2 <- FALSE
  
  while((iter < max.iter & !converged1) | !converged2){
    ## Update Beta
    # Find inverse of Sigma before
    Sigma_inv <- solve(Sigma)
    
    # Initialize product objects before loop
    Xt_SigInv_X <- matrix(0, ncol = length(Beta), nrow = length(Beta))
    Xt_SigInv_y <- rep(0, length(Beta))
    
    # Summing over observations
    for (i in 1:n) {
      Xt_SigInv_X <- Xt_SigInv_X + (t(X[[i]]) %*% Sigma_inv %*% X[[i]])
      Xt_SigInv_y <- Xt_SigInv_y + (t(X[[i]]) %*% Sigma_inv %*% y[[i]])
    }
    
    # Update
    BetaNew <- solve(Xt_SigInv_X) %*% Xt_SigInv_y
    
    ## Update Sigma
    # Initialize residual object
    residuals <- vector("list", n)
    
    # Sum over observations 
    for (i in 1:n) {
      residuals[[i]] <- y[[i]] - X[[i]] %*% BetaNew
    }
    
    # Update
    SigmaNew <- Reduce("+", lapply(residuals, function(r) r %*% t(r))) / n
    
    # Check for beta convergence
    if(norm(BetaNew - Beta, type = "2") < tol){
      converged1 <- TRUE
    }
    
    # Check for sigma convergence
    if(norm(SigmaNew - Sigma, type = "2") < tol){
      converged2 <- TRUE
    }
    
    # Keeping track of iterations
    iter <- iter + 1
    
    if(verbose == TRUE){
      cat("Iteration", iter, "completed \n")
    }
    
    # Copy for next iteration
    Beta <- BetaNew
    Sigma <- SigmaNew
  }
  
  # Check if maximum iterations were reached before convergence
  if (!converged1 | !converged2) {
    warning("Algorithm did not converge within maximum iterations.")
  }
  
  return(list(Beta = BetaNew, Sigma = SigmaNew))
}

## Initialize beta and Sigma
# Number of variables
p <- ncol(X[[1]])

# Initialize beta as a zero vector
beta_init <- rep(0, p)  

# Initialize Sigma as an identity matrix
Sigma_init <- diag(1, K)  

# Run it 
Algorithm(X, y, beta_init, Sigma_init, verbose = TRUE)

## Algorithm on data
# Reshaping Data
comp_case <- hw_data |>
  subset(
    !is.na(Week1) &
      !is.na(Week4) &
      !is.na(Week8)
  )
long_case <-
  stats::reshape(
    comp_case,
    direction = "long",
    varying = 3:5,
    sep = ""
  )[,-5]
names(long_case) <- c("Id","Treat","time","panss")

# Treat as factor and time as numeric (dunno)
long_case[["Treat"]] <- as.factor(long_case[["Treat"]])
# long_case$time <- as.numeric(long_case$time)

# Set up list of matrices for each observation
# Model matrices
X <- model.matrix(~ Treat * time, data = long_case)

# Split the rows of the model matrix by ID
X <- split(seq_len(nrow(long_case)), long_case[["Id"]]) |> 
  lapply(function(rows) X[rows, , drop = FALSE])

# Response vectors
y <- split(long_case, long_case[["Id"]]) |> 
  lapply(function(df) matrix(df[["panss"]], ncol = 1))

# Number of variables
p <- ncol(X[[1]])

# Initialize beta as a zero vector
beta_init <- rep(0, p)  

# Initialize Sigma as an identity matrix
Sigma_init <- diag(1, length(y[[1]]))  

# Run it
Output <- Algorithm(X, y, beta_init, Sigma_init, verbose = FALSE)
colnames(Output[["Sigma"]]) <- rownames(Output[["Sigma"]]) <- NULL
Output


### Question 2
library(compiler)
library(MCMCpack)
library(mvtnorm)
library(posterior)
library(parallel)

# Test Data
set.seed(123)
n <- 10000
p <- 5
K <- 3 
X <- list()
beta <- rnorm(p)
L <- matrix(rnorm(K * K), K, K)
Sigma <- L %*% t(L)
y <- list()
for (i in 1:n) {
  X[[i]] <- matrix(rnorm(p * K), K, p)
  y[[i]] <- X[[i]] %*% beta + MASS::mvrnorm(mu = rep(0, K), Sigma = Sigma)
}
data <- list(X = X, y = y)

MCMC.Fun <- function(data, p, K, 
                     mu0.vec = rep(0, p), sigma0.mat = diag(1, p), 
                     nu0 = 3, V0.mat = diag(1, K), 
                     B = 4, S = 500) {
  # Extract X and y from the input data
  X <- data$X
  y <- data$y
  n <- length(X)
  
  burn_in <- S / 2  # Burn-in period
  results <- list() # Store results
  
  for (chain in 1:B) {
    cat("Currently Running Chain", chain, "\n")
    
    # Initialize parameters
    Beta.vec <- rmvnorm(1, mean = mu0.vec, sigma = sigma0.mat)[1, ]
    Sigma.mat <- riwish(nu0, V0.mat)
    
    beta_samples <- matrix(NA, nrow = S, ncol = p)
    sigma_samples <- array(NA, dim = c(S, K, K))
    
    for (iter in 1:S) {
      # Pre-solve inverse matrices
      Inv.Sigma.mat <- solve(Sigma.mat)
      Inv.Sigma0.mat <- solve(sigma0.mat)
      
      # Step 1: Sample beta
      # Initialize beta_cov and beta_mean
      beta_cov <- matrix(0, nrow = p, ncol = p)
      beta_mean <- matrix(0, nrow = p, ncol = 1)
      
      for (j in 1:n) {
        beta_cov <- beta_cov + t(X[[j]]) %*% Inv.Sigma.mat %*% X[[j]]
        beta_mean <- beta_mean + t(X[[j]]) %*% Inv.Sigma.mat %*% y[[j]]
      }
      beta_cov <- solve(beta_cov + Inv.Sigma0.mat)
      beta_mean <- beta_cov %*% (beta_mean + Inv.Sigma0.mat %*% mu0.vec)
      
      # Sample Beta.vec
      Beta.vec <- as.vector(rmvnorm(1, mean = beta_mean, sigma = beta_cov))
      
      # # Step 1: Sample beta
      # beta_cov <- solve(Reduce('+', lapply(1:n, function(i) {
      #   t(X[[i]]) %*% solve(Sigma.mat) %*% X[[i]]
      # })) + solve(sigma0.mat))
      # 
      # beta_mean <- beta_cov %*% (Reduce('+', lapply(1:n, function(i) {
      #   t(X[[i]]) %*% solve(Sigma.mat) %*% y[[i]]
      # })) + solve(sigma0.mat) %*% mu0.vec)
      # 
      # Beta.vec <- as.vector(rmvnorm(1, mean = beta_mean, sigma = beta_cov))
      
      # Step 2: Sample Sigma
      residual_sum <- matrix(0, nrow = K, ncol = K) 
      for (i in 1:n) {
        residual <- y[[i]] - X[[i]] %*% Beta.vec 
        residual_sum <- residual_sum + residual %*% t(residual) 
      }
      
      Sigma.mat <- riwish(nu0 + n, V0.mat + residual_sum)
      
      # Store samples
      beta_samples[iter, ] <- Beta.vec
      sigma_samples[iter, , ] <- Sigma.mat
    }
    
    # Discard burn-in samples
    results[[chain]] <- list(
      beta = beta_samples[-(1:burn_in), ],
      sigma = sigma_samples[-(1:burn_in), , ]
    )
  }
  
  ChainNames <- paste("Chain", 1:B, sep = "")
  names(results) <- ChainNames
  
  # Reformatting for rhat() function use
  Beta.Array <- array(NA, dim = c(burn_in, B, p))
  for (i in 1:B) {
    for (j in 1:p) {
      Beta.Array[ , i, j] <- results[[ChainNames[i]]]$beta[, j]
    }
  }
  dimnames(Beta.Array) <- list(NULL,
                               ChainNames,
                               paste("Beta", 0:(p-1), sep = ""))
  
  # Compute Rhat statistics
  Beta.Rhats <- apply(Beta.Array, MARGIN = 3, FUN = rhat)
  
  # Return stuff
  return(list(
    results = results,
    rhats = Beta.Rhats
  ))
}

# Test function
result <- MCMC.Fun(data = data, p = p, K = K, S = 500, B = 4)

# Output Rhat statistics
result$rhats

# Parallel Function Version
MCMC.Parallel.Fun <- function(data, p, K, S = 500, B = 2, cores = detectCores(),
                              mu0.vec = rep(0, p), sigma0.mat = diag(1, p), 
                              nu0 = 3, V0.mat = diag(1, K)) {
  # Extract data
  X <- data$X
  y <- data$y
  n <- length(X)

  # Burn-in period
  burn_in <- S / 2  
  
  # Function to run a single chain
  run_chain <- function(chain_id) {
    cat("Running Chain", chain_id, "\n")
    
    Beta.vec <- rmvnorm(1, mean = mu0.vec, sigma = sigma0.mat)[1, ]
    Sigma.mat <- riwish(nu0, V0.mat)
    
    beta_samples <- matrix(NA, nrow = S, ncol = p)
    sigma_samples <- array(NA, dim = c(S, K, K))
    
    # Pre-solve inverse matrices
    Inv.Sigma0.mat <- solve(sigma0.mat)
    
    for (iter in 1:S) {
      # Pre-solve inverse matrices
      Inv.Sigma.mat <- solve(Sigma.mat)
      
      # Step 1: Sample beta
      # Initialize beta_cov and beta_mean
      beta_cov <- matrix(0, nrow = p, ncol = p)
      beta_mean <- matrix(0, nrow = p, ncol = 1)
      
      for (j in 1:n) {
        beta_cov <- beta_cov + t(X[[j]]) %*% Inv.Sigma.mat %*% X[[j]]
        beta_mean <- beta_mean + t(X[[j]]) %*% Inv.Sigma.mat %*% y[[j]]
      }
      beta_cov <- solve(beta_cov + Inv.Sigma0.mat)
      beta_mean <- beta_cov %*% (beta_mean + Inv.Sigma0.mat %*% mu0.vec)
      
      # Sample Beta.vec
      Beta.vec <- as.vector(rmvnorm(1, mean = beta_mean, sigma = beta_cov))
      
      # Step 2: Sample Sigma
      residual_sum <- matrix(0, nrow = K, ncol = K)
      for (i in 1:n) {
        residual <- y[[i]] - X[[i]] %*% Beta.vec
        residual_sum <- residual_sum + residual %*% t(residual)
      }
      Sigma.mat <- riwish(nu0 + n, V0.mat + residual_sum)
      
      # Store samples
      beta_samples[iter, ] <- Beta.vec
      sigma_samples[iter, , ] <- Sigma.mat
    }
    
    list(
      beta = beta_samples[-(1:burn_in), ],
      sigma = sigma_samples[-(1:burn_in), , ]
    )
  }
  
  # Run chains in parallel
  cl <- makeCluster(cores)  # Create a cluster with the specified number of cores
  clusterExport(cl, c("X", "y", "p", "K", "S", "mu0.vec", "sigma0.mat", 
                      "nu0", "V0.mat", "burn_in", "rmvnorm", "riwish"),
                envir = environment())
  results <- parLapply(cl, 1:B, run_chain)
  stopCluster(cl)  # Shut down the cluster
  
  # Reformatting for rhat() function
  ChainNames <- paste("Chain", 1:B, sep = "")
  names(results) <- ChainNames
  
  Beta.Array <- array(NA, dim = c(S - burn_in, B, p))
  for (i in 1:B) {
    for (j in 1:p) {
      Beta.Array[, i, j] <- results[[ChainNames[i]]]$beta[, j]
    }
  }
  dimnames(Beta.Array) <- list(NULL, ChainNames, paste("Beta", 0:(p - 1), sep = ""))
  
  # Univariate Rhat statistics
  Beta.Rhats <- apply(Beta.Array, MARGIN = 3, FUN = rhat)
  
  # Beta Quantiles
  Beta.Quantiles <- apply(Beta.Array, MARGIN = 3, FUN = quantile, prob = c(0.025, 0.975))
  
  # Return stuff
  return(list(
    results = results,
    rhats = Beta.Rhats,
    Beta.Quantiles = Beta.Quantiles
  ))
}

Comp.MCMC.Parallel.Fun <- cmpfun(MCMC.Parallel.Fun)

# Test function
outputs <- MCMC.Parallel.Fun(data, p = p, K = K, S = 500, B = detectCores(),
                             mu0.vec = rep(0, p), sigma0.mat = diag(1, p), 
                             nu0 = 3, V0.mat = diag(1, K))
outputs

# Format each unique element of covariance matrix
Pull.Cov <- function(gbg){
  Chain <- length(gbg$results)
  Dim <- dim(gbg$results$Chain1$sigma)[1]
  Ref.Mat <- diag(dim(gbg$results$Chain1$sigma)[2])
  RefCount <- sum(!upper.tri(Ref.Mat))
  Sigma.Array <- array(0, dim = c(Dim, RefCount, Chain))
  for(i in 1:Dim){
    for(j in 1:RefCount){
      for(k in 1:Chain)
        Sigma.Array[i, j, k] <- gbg$results[[paste("Chain", k, sep = "")]]$sigma[i,,][!upper.tri(Ref.Mat)][j]
    }
  }
  return(Sigma.Array)
}

# Univariate Rhat statistics for unique covariance elements
sigma.Rhats <- apply(Pull.Cov(outputs), MARGIN = 2, FUN = rhat)
sigma.Rhats

## Run in on dataset
# Set up list of matrices for each observation
# Model matrices
X <- model.matrix(~ Treat * time, data = long_case)

# Split the rows of the model matrix by ID
X <- split(seq_len(nrow(long_case)), long_case[["Id"]]) |> 
  lapply(function(rows) X[rows, , drop = FALSE])

# Response vectors
y <- split(long_case, long_case[["Id"]]) |> 
  lapply(function(df) matrix(df[["panss"]], ncol = 1))

# Number of variables
p <- ncol(X[[1]])
K <- nrow(X[[1]])

Data <- list(X = X, y = y)
# Run it
OutputS <- MCMC.Parallel.Fun(data = Data, p = p, K = K, S = 500, B = detectCores(),
                             mu0.vec = rep(0, p), sigma0.mat = diag(1, p), 
                             nu0 = 3, V0.mat = diag(1, K))
OutputS

# Univariate Rhat statistics for unique covariance elements
Sigma.Rhats <- apply(Pull.Cov(OutputS), MARGIN = 2, FUN = rhat)

Sigma.Rhats

# Creating correlation matrices
Corr.Array <- array(0, dim = c(250, 3, 3))
for(i in 1:250){
  Sig.Mat <- OutputS$results$Chain1$sigma[i, , ]
  D_inv_sqrt <- diag(1 / sqrt(diag(Sig.Mat)))
  Corr.Array[i, , ] <- D_inv_sqrt %*% Sig.Mat %*% D_inv_sqrt
}
Corr.Array

# Vector for each upper triangular rho
Rho.Mat <- matrix(0, nrow = 250, ncol = 3)

for(i in 1:250){
  # Rho_{1, 2}
  Rho.Mat[i, 1] <- Corr.Array[i, 2, 1]
  
  # Rho_{1, 3}
  Rho.Mat[i, 2] <- Corr.Array[i, 3, 1]
  
  # Rho_{2, 3}
  Rho.Mat[i, 3] <- Corr.Array[i, 3, 2]
}
colnames(Rho.Mat) <- c("Rho_{1, 2}", "Rho_{1, 3}", "Rho_{2, 3}")
Rho.Mat

# Finding 2.5% and 97.5% quantiles for each
quantile(Rho.Mat[, 1], prob = c(0.025, 0.975))
quantile(Rho.Mat[, 2], prob = c(0.025, 0.975))
quantile(Rho.Mat[, 3], prob = c(0.025, 0.975))
