### (ST 599) Homework 2
setwd("C:/Users/Ghcto/OneDrive/Desktop/School/Year 2/Winter 2025/ST 599/Homeworks")

# Load the data
library(Surrogate)
data("Schizo_PANSS")

# Subset the data to include only Week1, Week4, and Week8, including only complete cases, dropouts
# between week 4 and week 8, and dropouts between weeks 1 and 4:
hw_data <- Schizo_PANSS[,c("Id","Treat","Week1","Week4","Week8")] |>
  subset((!is.na(Week1) & !is.na(Week4) & !is.na(Week8)) |
           (!is.na(Week1) & !is.na(Week4) & is.na(Week8)) |
           (!is.na(Week1) & is.na(Week4) & is.na(Week8)))

# Function to find missingness patterns
gen_miss_patterns <- function(mat) {
  patterns <- mat |> is.na() |> apply(2, as.integer) |>
    apply(1, function(row) paste(row, collapse = "_"))
  return(patterns)
}

# Response Data
y_data <- hw_data[,c("Week1","Week4","Week8")]

# Find misingness patterns of the response data
miss_patterns <- gen_miss_patterns(y_data)

# Table and barplot of missingness patterns
tab_miss <- table(miss_patterns)
prop_miss <- prop.table(tab_miss)
prop_miss <- sort(prop_miss, decreasing = TRUE)
barplot(prop_miss, main = "Missingness patterns")

### Part (a)
## Simulate Data
set.seed(123)
n <- 1000
p <- 5
K <- 3
X <- list()
beta <- rnorm(p)
L <- matrix(rnorm(K * K), K, K)
Sigma <- L %*% t(L)
phi <- rnorm(p)
y <- list()
d <- rep(NA_integer_, n)
for (i in 1:n) {
  X[[i]] <- matrix(rnorm(p * K), K, p)
  y[[i]] <- X[[i]] %*% beta + MASS::mvrnorm(mu = rep(0, K), Sigma = Sigma)
  logit_p_m <- X[[i]] %*% phi + c(-2, -1, 1)
  p_m <- exp(logit_p_m) / sum(exp(logit_p_m))
  d_i <- rmultinom(1, 1, p_m)
  d[i] <- which(as.logical(d_i))
  y[[i]] <- y[[i]][1:d[i]]
}

## ECM Algorithm
# Conditional Expectation Function (Helper Function)
ConditionalExpectation <- function(y, X, Beta, Sigma, n) {
  # Pre-allocating space for moments
  FirstMoment <- vector("list", n)
  SecondMoment <- vector("list", n)
  Cond.Covar <- vector("list", n)
  
  # Calculating Expectations (row by row)
  for (i in seq_len(n)) {
    # Number of observed values & total possible values
    n0 <- length(y[[i]])         
    nk <- nrow(X[[i]])
    
    # No missing data, store observed values and outer product
    if (n0 == nk) {
      FirstMoment[[i]] <- y[[i]]
      SecondMoment[[i]] <- outer(y[[i]], y[[i]], FUN = "*")
      next
    }
    
    # Indices of observed and missing values
    obs_ind <- seq_len(n0)
    mis_ind <- setdiff(seq_len(nk), obs_ind)
    
    # Partition X based on observed/missing parts
    X_obs <- X[[i]][obs_ind, , drop = FALSE]
    X_mis <- X[[i]][mis_ind, , drop = FALSE]
    
    # Partition Sigma based on observed/missing parts
    Sigma_obs_obs <- Sigma[obs_ind, obs_ind, drop = FALSE] # Covariance among observed responses
    Sigma_mis_mis <- Sigma[mis_ind, mis_ind, drop = FALSE] # Covariance among missing responses
    Sigma_obs_mis <- Sigma[obs_ind, mis_ind, drop = FALSE] # Covariance between observed and missing responses
    Sigma_mis_obs <- Sigma[mis_ind, obs_ind, drop = FALSE] # Transpose of above
    
    # The expected values of the observed and missing responses
    mu_obs <- X_obs %*% Beta
    mu_mis <- X_mis %*% Beta
    
    # Compute E[Y_mis | y_obs]
    # Using multivariate normal conditional expectation formula:
    mu_mis_given_obs <- as.vector(mu_mis + Sigma_mis_obs %*% solve(Sigma_obs_obs) %*% (y[[i]] - mu_obs))
    
    # Compute Cov[Y_mis | y_obs]
    # Using multivariate normal conditional covariance formula:
    Cond.Cov <- Sigma_mis_mis - Sigma_mis_obs %*% solve(Sigma_obs_obs) %*% Sigma_obs_mis
    
    # First Moment Calculations
    # Full conditional mean vector
    # Filling in missing values with their conditional expectations
    FirMoment <- numeric(nk) # Initialize vector
    FirMoment[obs_ind] <- y[[i]]
    FirMoment[mis_ind] <- mu_mis_given_obs
    
    # Second Moment Calculations
    # Filling in missing values with their conditional expectations
    SecMoment <- matrix(0, nk, nk) # Initialize matrix
    # Outer product of observed values
    SecMoment[obs_ind, obs_ind] <- outer(y[[i]], y[[i]])            
    # Cross-product of observed and imputed missing values.
    SecMoment[obs_ind, mis_ind] <- outer(y[[i]], mu_mis_given_obs)  
    # Transpose of above
    SecMoment[mis_ind, obs_ind] <- outer(mu_mis_given_obs, y[[i]])  
    # Outer product of imputed missing values
    SecMoment[mis_ind, mis_ind] <- Cond.Cov + outer(mu_mis_given_obs, mu_mis_given_obs) 
    
    # Store the moments and conditional covariance
    FirstMoment[[i]] <- FirMoment
    SecondMoment[[i]] <- SecMoment
    Cond.Covar[[i]] <- Cond.Cov
  }
  
  return(list(FirstMoment = FirstMoment, SecondMoment = SecondMoment, Cond.Covar = Cond.Covar))
}

# Q Function (Helper Function)
Q.Fun <- function(y, X, Beta, Sigma, Beta.prev, Sigma.prev, Moments, n) {
  # Compute conditional moments using previous parameter estimates
  Moments <- ConditionalExpectation(y = y, X = X, Beta = Beta.prev, Sigma = Sigma.prev, n = n)
  
  # Pre-solve inverse
  Sigma_inv <- solve(Sigma)
  
  # Compute trace term
  trace_term <- 0
  
  for (i in seq_len(n)) {
    # Conditional Covariance Term
    Cov_yi <- Moments[["SecondMoment"]][[i]] - 
      Moments[["FirstMoment"]][[i]] %*% t(Moments[["FirstMoment"]][[i]])
    
    # Conditional expectation minus mu
    mean_diff <- Moments[["FirstMoment"]][[i]] - drop(X[[i]] %*% Beta)
    
    # Calculating the trace term
    trace_term <- trace_term + sum(diag((Cov_yi + mean_diff %*% t(mean_diff)) %*% Sigma_inv))
  }
  
  # Q
  Q <- -0.5 * n * log(det(Sigma)) - 0.5 * trace_term
  
  return(drop(Q))
}

# ECM Function
ECM.Fun <- function(X, y, Beta.Int, Sigma.Int, epsilon = 1e-7, max.iter = 50, verbose = FALSE){
  # Initialize Parameters
  Beta <- Beta.Int
  Sigma <- Sigma.Int
  n <- length(y)
  p <- ncol(X[[1]])
  k <- nrow(X[[1]])
  
  # Initialize function objects
  iter <- 0
  converged <- FALSE
  
  # Updating Steps
  while(iter < max.iter & !converged){
    ## E-Step, compute conditional expectations
    Moments <- ConditionalExpectation(y = y, X = X, Beta = Beta, Sigma = Sigma, n = n)

    ## Maximization Step
    # Update Beta
    # Pre-calc inverse of Sigma 
    SigmaInv <- solve(Sigma)
    
    # Initialize product objects before loop
    Obj1 <- matrix(0, ncol = p, nrow = p)
    Obj2 <- matrix(0, ncol = 1, nrow = p)
    
    # Summing over observations
    for (i in 1:n) {
      Obj1 <- Obj1 + (t(X[[i]]) %*% SigmaInv %*% X[[i]])
      Obj2 <- Obj2 + (t(X[[i]]) %*% SigmaInv %*% Moments[["FirstMoment"]][[i]])
    }
    
    # Update
    BetaNew <- solve(Obj1) %*% Obj2
    
    # Update Sigma
    # Using the residuals between the observed and predicted values weighted by their expected values
    
    # Initialize product objects before loop
    Obj3 <- matrix(0, ncol = k, nrow = k)
    
    for(i in 1:n){
      Obj3 <- Obj3 + Moments[["SecondMoment"]][[i]] - Moments[["FirstMoment"]][[i]] %*% t(X[[i]] %*% BetaNew) -
        t(Moments[["FirstMoment"]][[i]] %*% t(X[[i]] %*% BetaNew))  + (X[[i]] %*% BetaNew) %*% t(X[[i]] %*% BetaNew)
    }
    
    # Update
    SigmaNew <- Obj3 / n
    
    # Calculate Q Function Values
    Q1 <- Q.Fun(y = y, X = X, Beta = BetaNew, Sigma = SigmaNew, 
                Beta.prev = Beta, Sigma.prev = Sigma, 
                Moments = Moments, n = n)
    Q2 <- Q.Fun(y = y, X = X, Beta = Beta, Sigma = Sigma, 
                Beta.prev = Beta, Sigma.prev = Sigma, 
                Moments = Moments, n = n)
    
    # Check for convergence
    if(Q1 - Q2 < epsilon){
      converged <- TRUE
    }

    # Keeping track of iterations
    iter <- iter + 1

    # Printing
    if(verbose == TRUE){
      cat("Q New:", Q1, "Q Old:", Q2, "\n")
      cat("Iteration", iter, "completed \n")
    }

    # Copy for next iteration
    Beta <- BetaNew
    Sigma <- SigmaNew
  }

  # Check if maximum iterations were reached before convergence
  if (!converged) {
    warning("Algorithm did not converge within maximum iterations.")
  }
  
  # Output
  colnames(SigmaNew) <- NULL
  return(list(Beta = BetaNew, Sigma = SigmaNew, Iterations = iter))
}

# Initialize beta as a zero vector
beta_init <- rep(0, p)  

# Initialize Sigma as an identity matrix
Sigma_init <- diag(1, K)

## Run it
output <- ECM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                  epsilon = 1e-7, max.iter = 50, verbose = TRUE)
output

# Compute marginal asymptotic confidence intervals
MarginalCI <- function(y_list, X_list, SigmaMLE, BetaMLE, BetaTrue, alpha = 0.05, Verbose = TRUE) {
  # Calculate Moments
  n <- length(y_list)
  Moments <- ConditionalExpectation(y = y_list, X = X_list, Beta = BetaMLE, Sigma = SigmaMLE, n = n)
  
  # Pre-Compute Sigma Inverse
  InvSigma <- solve(SigmaMLE)
  
  # Compute the first summation
  sum1 <- 0
  for (Xi in X_list) {
    sum1 <- sum1 + t(Xi) %*% InvSigma %*% Xi
  }
  
  if(Verbose == TRUE){
    cat("First information term calculated \n")
    print(sum1)
  }
  
  # Compute the second summation
  sum2 <- 0
  for (i in seq_along(X_list)) {
    Xi <- X_list[[i]]
    Cov_yi <- Moments[["SecondMoment"]][[i]] - 
      Moments[["FirstMoment"]][[i]] %*% t(Moments[["FirstMoment"]][[i]])
    
    sum2 <- sum2 + t(Xi) %*% InvSigma %*% Cov_yi %*% InvSigma %*% Xi
  }
  
  if(Verbose == TRUE){
    cat("Second information term calculated \n")
    print(sum2)
  }
  
  # Compute C(y_i, X_i, β⋆, Σ⋆)
  C_list <- list()
  for (i in seq_along(X_list)) {
    mean_term <- Moments[["FirstMoment"]][[i]] - X_list[[i]] %*% BetaMLE
    Cov_yi <- Moments[["SecondMoment"]][[i]] - 
      Moments[["FirstMoment"]][[i]] %*% t(Moments[["FirstMoment"]][[i]])
    C_list[[i]] <- Cov_yi + mean_term %*% t(mean_term)
  }
  
  # Compute the third summation
  sum3 <- 0
  for (i in seq_along(X_list)) {
    Xi <- X_list[[i]]
    sum3 <- sum3 + t(Xi) %*% InvSigma %*% C_list[[i]] %*% InvSigma %*% Xi
  }
  
  if(Verbose == TRUE){
    cat("Third information term calculated \n")
    print(sum3)
  }
  
  # Compute the fourth summation
  sum4 <- 0
  for (i in seq_along(X_list)) {
    for (j in seq_along(X_list)) {
      if (i != j) {
        Xi <- X_list[[i]]
        Xj <- X_list[[j]]
        mu_i <- Moments[["FirstMoment"]][[i]] - Xi %*% BetaMLE
        mu_j <- Moments[["FirstMoment"]][[j]] - Xj %*% BetaMLE
        sum4 <- sum4 + t(Xi) %*% InvSigma %*% (mu_i %*% t(mu_j)) %*% InvSigma %*% Xj
      }
    }
  }
  
  if(Verbose == TRUE){
    cat("Forth information term calculated \n")
    print(sum4)
  }
  
  # Final result for information matrix
  Info <- sum1 - sum2 - sum3 - sum4
  
  # Compute Confidence Intervals
  CI.Mat <- matrix(NA, nrow = length(BetaTrue), ncol = 2)
  SE.Vec <- rep(0, length(BetaTrue))
  Zcrit <- qnorm(1 - alpha/2)
  InvInfo <- solve(Info)
  for(i in 1:length(BetaTrue)){
    SE.Vec[i] <- sqrt(InvInfo[i, i])
    CI.Mat[i, ] <- BetaMLE[i] + c(-1, 1)*SE.Vec[i]
  }
  
  # Checking if MLE is in the confidence interval
  Cap <- (BetaTrue >= CI.Mat[, 1]) & (BetaTrue <= CI.Mat[, 2])
  
  CI <- cbind.data.frame(CI.Mat, Cap)
  rownames(CI) <- paste("Beta", seq_along(BetaTrue))
  colnames(CI) <- c("Lower", "Upper", "Captured")
  
  return(list(Information = Info,
              ConfidenceIntervals = CI,
              StandardErrors = SE.Vec))
}

# Run it
MarginalCI(y_list = y, X_list = X, 
           SigmaMLE = output$Sigma, BetaMLE = output$Beta, BetaTrue = beta,
           alpha = 0.05)

### Part (b)
# Turn each row into a list element 
Y_Data <- split(y_data, seq(nrow(y_data)))

# Remove NAs from each list element
Y_Data <- lapply(Y_Data , function(x) x[!is.na(x)])

# Reshaping Data
long_case <-
  stats::reshape(
    hw_data,
    direction = "long",
    varying = 3:5,
    sep = ""
  )[,-5]
names(long_case) <- c("Id","Treat","time","panss")

# Treat as factor
long_case[["Treat"]] <- as.factor(long_case[["Treat"]])

# Set up list of matrices for each observation
X.Mat <- model.matrix(~ Treat * time, data = long_case)

# Split the rows of the model matrix by ID
X.Mat <- split(seq_len(nrow(long_case)), long_case[["Id"]]) |> 
  lapply(function(rows) X.Mat[rows, , drop = FALSE])

# Initialize as a zero vector
beta_init <- rep(0, ncol(X.Mat[[1]]))  

# Initialize Sigma as an identity matrix
Sigma_init <- diag(1, nrow(X.Mat[[1]]))

## Run it
OutputData <- ECM.Fun(X = X.Mat, y = Y_Data, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                   epsilon = 1e-7, max.iter = 50, verbose = TRUE)
OutputData

# The captured elements are useless but CI and SE can be extracted
MarginalCI(y_list = Y_Data, X_list = X.Mat, 
           SigmaMLE = OutputData$Sigma, BetaMLE = OutputData$Beta, BetaTrue = beta[-1],
           alpha = 0.05)

### Part (c)
EM.Fun <- function(X, y, Beta.Int, Sigma.Int, epsilon = 1e-7, max.iter = 100, verbose = FALSE){
  # Initialize Parameters
  Beta <- Beta.Int
  Sigma <- Sigma.Int
  n <- length(y)
  p <- ncol(X[[1]])
  k <- nrow(X[[1]])
  
  # Initialize function objects
  iter <- 0
  Iter <- 0
  converged3 <- FALSE
  
  # Updating Steps
  while(iter < max.iter & !converged3){
    # Initialize function objects
    converged2 <- converged1 <- FALSE
    
    ## E-Step, compute conditional expectations
    Moments <- ConditionalExpectation(y = y, X = X, Beta = Beta, Sigma = Sigma, n = n)
    
    ## Maximization Step
    BetaS <- Beta
    SigmaS <- Sigma
    # Iterate until convergence at each iteration
    while((Iter < max.iter & !converged1) | (Iter < max.iter & !converged2)){
      # Update Beta
      # Pre-calc inverse of Sigma 
      SigmaInv <- solve(SigmaS)
      
      # Initialize product objects before loop
      Obj1 <- matrix(0, ncol = p, nrow = p)
      Obj2 <- matrix(0, ncol = 1, nrow = p)
      
      # Summing over observations
      for (i in 1:n) {
        Obj1 <- Obj1 + (t(X[[i]]) %*% SigmaInv %*% X[[i]])
        Obj2 <- Obj2 + (t(X[[i]]) %*% SigmaInv %*% Moments[["FirstMoment"]][[i]])
      }
      
      # Update
      BetaSNew <- solve(Obj1) %*% Obj2
      
      # Update Sigma
      # Using the residuals between the observed and predicted values weighted by their expected values
      # Initialize product objects before loop
      Obj3 <- matrix(0, ncol = k, nrow = k)
      for(i in 1:n){
        Obj3 <- Obj3 + Moments[["SecondMoment"]][[i]] - Moments[["FirstMoment"]][[i]] %*% t(X[[i]] %*% BetaSNew) -
          t(Moments[["FirstMoment"]][[i]] %*% t(X[[i]] %*% BetaSNew)) + (X[[i]] %*% BetaSNew) %*% t(X[[i]] %*% BetaSNew)
      }
      
      # Update
      SigmaSNew <- Obj3 / n
      
      # Check for beta convergence
      if(norm(BetaSNew - BetaS, type = "2") < epsilon){
        converged1 <- TRUE
      }
      
      # Check for sigma convergence
      if(norm(SigmaSNew - SigmaS, type = "2") < epsilon){
        converged2 <- TRUE
      }
      
      # Keeping track of iterations
      Iter <- Iter + 1
      
      # Printing
      if(verbose == TRUE){
        cat("BetaSNew", BetaSNew, "\n")
        cat("SigmaSNew:", "\n")
        print(SigmaSNew)
        cat("S Iteration", Iter, "completed \n")
      }
      
      # Copy for next iteration
      BetaS <- BetaSNew
      SigmaS <- SigmaSNew
    }
    
    # Copy parameters over
    BetaNew <- BetaS
    SigmaNew <- SigmaS
    
    # Calculate Q Function Values
    Q1 <- Q.Fun(y = y, X = X, Beta = BetaNew, Sigma = SigmaNew, 
                Beta.prev = Beta, Sigma.prev = Sigma, 
                Moments = Moments, n = n)
    Q2 <- Q.Fun(y = y, X = X, Beta = Beta, Sigma = Sigma, 
                Beta.prev = Beta, Sigma.prev = Sigma, 
                Moments = Moments, n = n)
    
    # Check for convergence
    if(Q1 - Q2 < epsilon){
      converged3 <- TRUE
    }
    
    # Keeping track of iterations
    iter <- iter + 1
    
    # Printing
    if(verbose == TRUE){
      cat("Q New:", Q1, "Q Old:", Q2, "\n")
      cat("T Iteration", iter, "completed \n")
    }
    
    # Copy for next iteration
    Beta <- BetaNew
    Sigma <- SigmaNew
  }
  
  # Check if maximum iterations were reached before convergence
  if (!converged3) {
    warning("Algorithm did not converge within maximum iterations.")
  }
  
  # Output
  colnames(SigmaNew) <- NULL
  return(list(Beta = BetaNew, Sigma = SigmaNew,
              Iterations_S = Iter, Iterations_T = iter))
}

# Initialize beta as a zero vector
beta_init <- rep(0, p)  

# Initialize Sigma as an identity matrix
Sigma_init <- diag(1, K)  

## Run it
output2 <- EM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                  epsilon = 1e-7, max.iter = 50, verbose = TRUE)
output2

# Compute marginal asymptotic confidence intervals
MarginalCI(y_list = y, X_list = X, 
           SigmaMLE = output2$Sigma, BetaMLE = output2$Beta, BetaTrue = beta)


# Parallel Package
library(future.apply)
library(progressr)
handlers("txtprogressbar")  # Enables a visible progress bar

# Function to assess coverage
Coverage.Fun <- function(seedint, S, alpha, BetaTrue, phi, K, n = 1000) {
  # Number of parameters
  p <- length(BetaTrue)
  
  # Progress bar
  with_progress({
    prog <- progressor(along = 1:S)
    
    # Define function to run a single iteration
    run_iteration <- function(i) {
      
      set.seed(seedint + i) # Set new seed
      
      X <- vector("list", n)
      y <- vector("list", n)
      d <- integer(n)
      
      for (j in 1:n) {
        X[[j]] <- matrix(rnorm(p * K), K, p)
        y[[j]] <- X[[j]] %*% BetaTrue + MASS::mvrnorm(mu = rep(0, K), Sigma = Sigma)
        logit_p_m <- X[[j]] %*% phi + c(-2, -1, 1)
        p_m <- exp(logit_p_m) / sum(exp(logit_p_m))
        d_j <- rmultinom(1, 1, p_m)
        d[j] <- which(as.logical(d_j))
        y[[j]] <- y[[j]][1:d[j]]
      }
      
      beta_init <- rep(0, p)
      Sigma_init <- diag(K)
      
      # ECM
      start_time_ecm <- Sys.time()
      Output <- ECM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                        epsilon = 1e-4, max.iter = 50, verbose = FALSE)
      end_time_ecm <- Sys.time()
      ecm_runtime <- as.numeric(difftime(end_time_ecm, start_time_ecm, units = "secs"))
      
      # EM
      start_time_em <- Sys.time()
      Output2 <- EM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                        epsilon = 1e-4, max.iter = 50, verbose = FALSE)
      end_time_em <- Sys.time()
      em_runtime <- as.numeric(difftime(end_time_em, start_time_em, units = "secs"))
      
      # Marginal CI
      ECMCI <- MarginalCI(y_list = y, X_list = X, 
                          SigmaMLE = Output$Sigma, BetaMLE = Output$Beta, BetaTrue = BetaTrue,
                          alpha = alpha, Verbose = FALSE)
      EMCI <- MarginalCI(y_list = y, X_list = X, 
                         SigmaMLE = Output2$Sigma, BetaMLE = Output2$Beta, BetaTrue = BetaTrue,
                         alpha = alpha, Verbose = FALSE)
      
      prog()  # Update the progress bar
      
      return(list(
        ECM_Coverage = ECMCI$ConfidenceIntervals$Captured,
        EM_Coverage = EMCI$ConfidenceIntervals$Captured,
        ECM_Speed = c(ecm_runtime, Output$Iterations, 0),
        EM_Speed = c(em_runtime, Output2$Iterations_T, Output2$Iterations_S)
      ))
      
    }
    
    # Set up parallel backend
    plan(multisession, workers = availableCores())
    
    # Run iterations in parallel
    results <- future_lapply(1:S, run_iteration, future.seed = TRUE)
  })
  
  # Extract and combine results
  Coverage.Mat1 <- do.call(rbind, lapply(results, `[[`, "ECM_Coverage"))
  Coverage.Mat2 <- do.call(rbind, lapply(results, `[[`, "EM_Coverage"))
  TimeIterMat1 <- do.call(rbind, lapply(results, `[[`, "ECM_Speed"))
  TimeIterMat2 <- do.call(rbind, lapply(results, `[[`, "EM_Speed"))
  
  # Compute final statistics
  Speed1 <- colMeans(TimeIterMat1)
  Speed2 <- colMeans(TimeIterMat2)
  names(Speed1) <- names(Speed2) <- c("Avg_Run_Sec", "Avg_t_iter", "Avg_s_iter")
  
  Coverage1 <- colMeans(Coverage.Mat1)
  Coverage2 <- colMeans(Coverage.Mat2)
  names(Coverage1) <- names(Coverage2) <- paste("Beta", seq_along(BetaTrue))
  
  return(list(ECM_Coverage = Coverage1,
              EM_Coverage = Coverage2,
              ECM_Speed = Speed1,
              EM_Speed = Speed2))
}

# Run it
Coverage.Fun(seedint = 123, S = 500, alpha = 0.05, BetaTrue = beta, phi = phi, K = K, n = 5000)

# Without parallelization
Coverage.FunOld <- function(seedint, S, alpha, BetaTrue, phi, K, n = 1000){
  Coverage.Mat2 <- Coverage.Mat1 <- matrix(NA, nrow = S, ncol = length(BetaTrue))
  TimeIterMat2 <- TimeIterMat1 <- matrix(NA, nrow = S, ncol = 3)
  for(i in 1:S){
    # Changing the seed each iteration
    set.seed(seedint + i)
    
    # Simulate Data
    p <- length(BetaTrue)
    X <- list()
    y <- list()
    d <- rep(NA_integer_, n)
    for (j in 1:n) {
      X[[j]] <- matrix(rnorm(p * K), K, p)
      y[[j]] <- X[[j]] %*% BetaTrue + MASS::mvrnorm(mu = rep(0, K), Sigma = Sigma)
      logit_p_m <- X[[j]] %*% phi + c(-2, -1, 1)
      p_m <- exp(logit_p_m) / sum(exp(logit_p_m))
      d_j <- rmultinom(1, 1, p_m)
      d[j] <- which(as.logical(d_j))
      y[[j]] <- y[[j]][1:d[j]]
    }
    
    # Initialize parameters for algorithm
    beta_init <- rep(0, p)
    Sigma_init <- diag(K)
    
    # Perform ECM
    start_time_ecm <- Sys.time()
    Output <- ECM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                      epsilon = 1e-4, max.iter = 50, verbose = FALSE)
    end_time_ecm <- Sys.time()
    ecm_runtime <- as.numeric(difftime(end_time_ecm, start_time_ecm, units = "secs"))
    
    # Perform EM
    start_time_em <- Sys.time()
    Output2 <- EM.Fun(X, y, Beta.Int = beta_init, Sigma.Int = Sigma_init, 
                      epsilon = 1e-4, max.iter = 50, verbose = FALSE)
    end_time_em <- Sys.time()
    em_runtime <- as.numeric(difftime(end_time_em, start_time_em, units = "secs"))
    
    # Calculate the marginal asymptotic confidence intervals for both
    ECMCI <- MarginalCI(y_list = y, X_list = X, 
                        SigmaMLE = Output$Sigma, BetaMLE = Output$Beta, BetaTrue = BetaTrue,
                        alpha = alpha, Verbose = FALSE)
    EMCI <- MarginalCI(y_list = y, X_list = X, 
                       SigmaMLE = Output2$Sigma, BetaMLE = Output2$Beta, BetaTrue = BetaTrue,
                       alpha = alpha, Verbose = FALSE)
    
    # Save runtime and iteration results
    # ECM
    TimeIterMat1[i, 1] <- ecm_runtime
    TimeIterMat1[i, 2] <- Output$Iterations
    TimeIterMat1[i, 3] <- 0
    
    # EM
    TimeIterMat2[i, 1] <- em_runtime
    TimeIterMat2[i, 2] <- Output2$Iterations_T
    TimeIterMat2[i, 3] <- Output2$Iterations_S
    
    # Save capture results
    Coverage.Mat1[i, ] <- ECMCI$ConfidenceIntervals$Captured
    Coverage.Mat2[i, ] <- EMCI$ConfidenceIntervals$Captured
    print(i)
  }
  # Obtain average Runtime and iterations for both methods
  Speed1 <- apply(TimeIterMat1, 2, FUN = mean)
  Speed2 <- apply(TimeIterMat2, 2, FUN = mean)
  names(Speed2) <- names(Speed1) <- c("Avg_Run_Sec", "Avg_t_iter", "Avg_s_iter")
  
  # Calculate the proportion of intervals that contain the true beta
  Coverage1 <- apply(Coverage.Mat1, 2, FUN = mean)
  Coverage2 <- apply(Coverage.Mat2, 2, FUN = mean)
  names(Coverage2) <- names(Coverage1) <- paste("Beta", seq_along(BetaTrue))
  return(list(ECM_Coverage = Coverage1,
              EM_Coverage = Coverage2,
              ECM_Speed = Speed1,
              EM_Speed = Speed2))
}
Coverage.FunOld(seedint = 123, S = 10, alpha = 0.05, BetaTrue = beta, phi = phi, K = K)