### Functions and packages
# Load packages
packages <- c("FNN", "glmnet", "leaps", "future.apply","Rcpp", "RcppArmadillo", "Amelia")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Cpp Functions
sourceCpp("./Missing Data Project/C Code/BIMI.cpp")
# sourceCpp("BIMI.cpp")

# Set up progress bar :)
# library(progressr)
# handlers(global = TRUE)
# handlers("txtprogressbar")

### Generate Data
GenData.Fun <- function(n, Beta, Sigma, ErrorVar, 
                        MissMech = "MCAR", MissProp = 0.20, Resp_Miss = TRUE) {
  # Number of predictors
  p <- length(Beta) - 1
  
  # Generate predictor matrix X with specified covariance structure Sigma
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Generate response variables y using X and Beta
  epsilon <- rnorm(n, mean = 0, sd = sqrt(ErrorVar))
  y <- Beta[1] + X %*% Beta[-1] + epsilon
  
  # Remove observations based on type of missingness mechanism
  if (MissMech == "MCAR") {
    # MCAR: Missingness is independent of the data
    miss_ind <- matrix(runif(n * p) < MissProp, n, p)
    X[miss_ind] <- NA
    
    # If response is missing
    if (Resp_Miss == TRUE){
      y[runif(n) < MissProp] <- NA
    }
    
  } else if (MissMech == "MAR") {
    # If response is missing
    if (Resp_Miss == TRUE){
      # MAR: Missingness depends on a subset of variables
      Dataset <- cbind(X = X, y = y)
      
      # Initialize missingness matrix
      Miss.Mat <- matrix(0, nrow = dim(Dataset)[1], ncol = dim(Dataset)[2])
      
      # What are the observed columns
      NonMissVars <- Dataset[, (round((p + 1) / 2) + 1):p] 
      
      # PCA
      NonMissVars.pc <- princomp(NonMissVars, scores = TRUE)$scores[, 1]
      
      for (j in 1:round((p + 1) / 2)) {
        # Scaled to account for half of the columns being fully observed
        MissProbVec <- 4 * MissProp * pnorm(NonMissVars.pc)
        
        # Indices 
        MissInd <- rbinom(n, 1, MissProbVec)
        
        Dataset[MissInd == 1, j] <- NA
      }
      # Put y column back into first slot
      Dataset <- Dataset[, c(ncol(Dataset), 1:(ncol(Dataset) - 1))]
      
    } else {
      # MAR: Missingness depends on a subset of variables
      # Initialize missingness matrix
      Miss.Mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
      
      # What are the observed columns
      NonMissVars <- X[, (round((p + 1) / 2) + 1):p] 
      
      # PCA
      NonMissVars.pc <- princomp(NonMissVars, scores = TRUE)$scores[, 1]
      
      for (j in 1:round((p + 1) / 2)) {
        # Scaled to account for half of the columns being fully observed
        MissProbVec <- 4 * MissProp * pnorm(NonMissVars.pc)
        
        # Indices 
        MissInd <- rbinom(n, 1, MissProbVec)
        
        X[MissInd == 1, j] <- NA
      }
    }
    
  } else if (MissMech == "MNAR") {
    # If response is missing
    if (Resp_Miss == TRUE){
      Dataset <- cbind(y = y, X = X)
      # MNAR: Missingness depends on the magnitude of the observed values
      for (j in 1:(p+1)) {
        # Use the normal CDF to calculate probability of missingness based on the magnitude of variable
        miss_probs <- 2 * MissProp * pnorm(Dataset[, j], mean = mean(Dataset[, j]), sd = sd(Dataset[, j]))
        
        # Apply missingness based on the calculated probabilities
        Dataset[, j][runif(n) < miss_probs] <- NA
      }
      
    } else {
      # MNAR: Missingness depends on the magnitude of the observed values
      for (j in 1:p) {
        # Use the normal CDF to calculate probability of missingness based on the magnitude of X
        miss_probs <- 2 * MissProp * pnorm(X[, j], mean = mean(X[, j]), sd = sd(X[, j]))
        
        # Apply missingness based on the calculated probabilities
        X[, j][runif(n) < miss_probs] <- NA
      }
    }
    
  } else {
    stop("Missingness mechanism not supported. Choose either 'MCAR', 'MAR', or 'MNAR'.")
  }
  
  if (Resp_Miss == TRUE) {
    if(MissMech == "MCAR"){
      # Combine
      Dataset <- cbind(y = y, X = X)
    } 
    
    # Remove any rows with all X's NAs
    NACheck <- apply(Dataset, 1, function(x) any(!is.na(x)))
    Dataset <- Dataset[NACheck, ]
    
  } else {
    # Remove any rows with all X's NAs
    NACheck <- apply(X, 1, function(x) any(!is.na(x)))
    X <- X[NACheck, ]
    y <- y[NACheck, ]
    
    # Return the dataset
    Dataset <- cbind(y = y, X = X)
  }
  
  colnames(Dataset) <- c("Y", paste("X", 1:(length(Beta)-1), sep = ""))
  
  return(Dataset)
}

data_generator <- function(n, Sigma, marginal_dists, parameters) {
  # Number of predictors and basic check
  p <- length(marginal_dists)
  if (length(parameters) != p) stop("Length of parameters must equal length of marginal disrtibutions provided")
  
  # Step 1: Generate multivariate normal latent variables
  Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Step 2: Transform each variable to the desired marginal
  X <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    # Percentiles under standard normal
    p_score <- pnorm(Z[, j])
    
    # Build the quantile function call
    dist_name <- marginal_dists[[j]]     # e.g., "norm", "exp", "t"
    qfun_name <- paste0("q", dist_name)
    
    if (!exists(qfun_name)) stop(paste("Quantile function", qfun_name, "does not exist"))
    
    qfun <- get(qfun_name)
    
    # Call quantile function with parameters
    X[, j] <- do.call(qfun, c(list(p = p_score), parameters[[j]]))
  }
  
  colnames(X) <- paste0("X", 1:p)
  return(X)
}

# Combination of the two above functions
OldDataGeneration <- function(n, Beta, Sigma, ErrorVar, error_dist = NULL, error_parms = NULL,
                            marginal_dists = NULL, parameters = NULL,
                            MissMech = "MCAR", MissProp = 0.20, Resp_Miss = TRUE) {
  ### Data Generation
  if (is.null(marginal_dists) || is.null(parameters)) {
    # Number of predictors
    p <- length(Beta) - 1
    
    # Generate predictor matrix X with specified covariance structure Sigma
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
  } else {
    # Number of predictors and basic check
    p <- length(marginal_dists)
    if (length(parameters) != p) stop("Length of parameters must equal length of marginal disrtibutions provided")
    
    # Step 1: Generate multivariate normal latent variables
    Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    # Step 2: Transform each variable to the desired marginal
    X <- matrix(NA, nrow = n, ncol = p)
    for (j in 1:p) {
      # Percentiles under standard normal
      p_score <- pnorm(Z[, j])
      
      # Build the quantile function call
      dist_name <- marginal_dists[[j]]     # e.g., "norm", "exp", "t"
      qfun_name <- paste0("q", dist_name)
      
      if (!exists(qfun_name)) stop(paste("Quantile function", qfun_name, "does not exist"))
      
      qfun <- get(qfun_name)
      
      # Call quantile function with parameters
      X[, j] <- do.call(qfun, c(list(p = p_score), parameters[[j]]))
    }
  }
  
  # Generate error terms
  if (is.null(error_dist) || is.null(error_parms)) {
    epsilon <- rnorm(n, mean = 0, sd = sqrt(ErrorVar))
  } else {
    error_fun <- paste0("r", error_dist)
    
    if (!exists(error_fun)) stop(paste("Random data generating function", error_fun, "does not exist"))
    
    error_fun <- get(paste0("r", error_dist))
    
    epsilon <- do.call(error_fun, c(list(n = n), error_parms))
  }
  
  # Generate response variables y using X and Beta
  y <- Beta[1] + X %*% Beta[-1] + epsilon
  
  ## No Missingness
  if(MissMech == "None") {
    # Combine
    Dataset <- cbind(y = y, X = X)
    
    colnames(Dataset) <- c("Y", paste("X", 1:(length(Beta) - 1), sep = ""))
    
    return(Dataset)
  }
  
  ### Missingness
  # Remove observations based on type of missingness mechanism
  if (MissMech == "MCAR") {
    # MCAR: Missingness is independent of the data
    miss_ind <- matrix(runif(n * p) < MissProp, n, p)
    X[miss_ind] <- NA
    
    # If response is missing
    if (Resp_Miss == TRUE){
      y[runif(n) < MissProp] <- NA
    }
    
  } else if (MissMech == "MAR") {
    # If response is missing
    if (Resp_Miss == TRUE){
      # MAR: Missingness depends on a subset of variables
      Dataset <- cbind(X = X, y = y)
      
      # Initialize missingness matrix
      Miss.Mat <- matrix(0, nrow = dim(Dataset)[1], ncol = dim(Dataset)[2])
      
      # What are the observed columns
      NonMissVars <- Dataset[, (round((p + 1) / 2) + 1):p] 
      
      # PCA
      NonMissVars.pc <- princomp(NonMissVars, scores = TRUE)$scores[, 1]
      
      for (j in 1:round((p + 1) / 2)) {
        # Scaled to account for half of the columns being fully observed
        MissProbVec <- 4 * MissProp * pnorm(NonMissVars.pc)
        
        # Indices 
        MissInd <- rbinom(n, 1, MissProbVec)
        
        Dataset[MissInd == 1, j] <- NA
      }
      # Put y column back into first slot
      Dataset <- Dataset[, c(ncol(Dataset), 1:(ncol(Dataset) - 1))]
      
    } else {
      # MAR: Missingness depends on a subset of variables
      # Initialize missingness matrix
      Miss.Mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
      
      # What are the observed columns
      NonMissVars <- X[, (round((p + 1) / 2) + 1):p] 
      
      # PCA
      NonMissVars.pc <- princomp(NonMissVars, scores = TRUE)$scores[, 1]
      
      for (j in 1:round((p + 1) / 2)) {
        # Scaled to account for half of the columns being fully observed
        MissProbVec <- 4 * MissProp * pnorm(NonMissVars.pc)
        
        # Indices 
        MissInd <- rbinom(n, 1, MissProbVec)
        
        X[MissInd == 1, j] <- NA
      }
    }
    
  } else if (MissMech == "MNAR") {
    # If response is missing
    if (Resp_Miss == TRUE){
      Dataset <- cbind(y = y, X = X)
      # MNAR: Missingness depends on the magnitude of the observed values
      for (j in 1:(p+1)) {
        # Use the normal CDF to calculate probability of missingness based on the magnitude of variable
        miss_probs <- 2 * MissProp * pnorm(Dataset[, j], mean = mean(Dataset[, j]), sd = sd(Dataset[, j]))
        
        # Apply missingness based on the calculated probabilities
        Dataset[, j][runif(n) < miss_probs] <- NA
      }
      
    } else {
      # MNAR: Missingness depends on the magnitude of the observed values
      for (j in 1:p) {
        # Use the normal CDF to calculate probability of missingness based on the magnitude of X
        miss_probs <- 2 * MissProp * pnorm(X[, j], mean = mean(X[, j]), sd = sd(X[, j]))
        
        # Apply missingness based on the calculated probabilities
        X[, j][runif(n) < miss_probs] <- NA
      }
    }
    
  } else {
    stop("Missingness mechanism not supported. Choose either 'MCAR', 'MAR', or 'MNAR'.")
  }
  
  if (Resp_Miss == TRUE) {
    if(MissMech == "MCAR"){
      # Combine
      Dataset <- cbind(y = y, X = X)
    } 
    
    # Remove any rows with all X's NAs
    NACheck <- apply(Dataset, 1, function(x) any(!is.na(x)))
    Dataset <- Dataset[NACheck, ]
    
  } else {
    # Remove any rows with all X's NAs
    NACheck <- apply(X, 1, function(x) any(!is.na(x)))
    X <- X[NACheck, ]
    y <- y[NACheck, ]
    
    # Return the dataset
    Dataset <- cbind(y = y, X = X)
  }
  
  colnames(Dataset) <- c("Y", paste("X", 1:(length(Beta)-1), sep = ""))
  
  return(Dataset)
}

#' Current data generation function to properly implement missingness
#'
#' @param n Sample size
#' @param Beta Coefficient vector
#' @param Sigma Covariance matrix
#' @param ErrorVar Residual error variance
#' @param error_dist Residual error distribution 
#' @param error_parms Residual error distribution parameters
#' @param marginal_dists Covariate distributions
#' @param parameters # Covariate distributions' parameters
#' @param marginal_means # Covariate population means
#' @param MissProp Target missingness proportions (1 - proportion of complete cases)
#' @param Eta  Weighting vector (without intercept if not found) determining the missingness mechanism, NULL or "None" for no missingness
#' @param MissCol Variable / column to induce missingness in eg. "X1" or 2
#'
#' @returns
#' @export
#'
#' @examples
DataEtaBetaZeta <- function(n, Beta, Sigma, ErrorVar, error_dist = NULL, error_parms = NULL,
                            marginal_dists = NULL, parameters = NULL, marginal_means = NULL,
                            MissProp = 0.20, Eta = NULL, MissCol, MissMechType = "logit") {
  ### Data Generation
  if (is.null(marginal_dists) || is.null(parameters)) {
    print("No covariate distributions and or parmeters given. Normal covariates will be generated.")
    
    # Number of predictors
    p <- length(Beta) - 1
    
    # Generate predictor matrix X with specified covariance structure Sigma
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
  } 
  else {
    # Number of predictors and basic check
    p <- length(marginal_dists)
    if (length(parameters) != p) stop("Length of parameters must equal length of marginal disrtibutions provided")
    
    # Step 1: Generate multivariate normal latent variables
    Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    # Step 2: Transform each variable to the desired marginal
    X <- matrix(NA, nrow = n, ncol = p)
    for (j in 1:p) {
      # Percentiles under standard normal
      p_score <- pnorm(Z[, j])
      
      # Build the quantile function call
      dist_name <- marginal_dists[[j]]     # e.g., "norm", "exp", "t"
      qfun_name <- paste0("q", dist_name)
      
      if (!exists(qfun_name)) stop(paste("Quantile function", qfun_name, "does not exist"))
      
      qfun <- get(qfun_name)
      
      # Call quantile function with parameters
      X[, j] <- do.call(qfun, c(list(p = p_score), parameters[[j]]))
    }
  }
  
  # Generate error terms
  if (is.null(error_dist) || is.null(error_parms)) {
    print("No error distribution and or parmeters given. Normal errors will be generated.")
    
    epsilon <- rnorm(n, mean = 0, sd = sqrt(ErrorVar))
  } 
  else {
    error_fun <- paste0("r", error_dist)
    
    if (!exists(error_fun)) stop(paste("Random data generating function", error_fun, "does not exist"))
    
    error_fun <- get(paste0("r", error_dist))
    
    epsilon <- do.call(error_fun, c(list(n = n), error_parms))
  }
  
  # Generate response variables y using X and Beta
  y <- Beta[1] + X %*% Beta[-1] + epsilon
  
  ## No Missingness
  if (is.null(Eta) | any(Eta == "None") | MissMechType == "None" | MissProp == 0) {
    # Combine
    Dataset <- cbind(y = y, X = X)
    
    colnames(Dataset) <- c("Y", paste("X", 1:(length(Beta) - 1), sep = ""))

    return(Dataset)
  }
  ### Missingness
  if (MissMechType == "logit") {
    # If intercept is given, no need to find it
    if (length(Eta) != p + 2) {
      # Grid search for intercept (Eta[1]) to obtain target missingness proportion
      Eta <- c(0, Eta)
      
      popMat <- OldDataGeneration(n = 1e6, # 1e8 amount blew up my PC
                                  Beta = Beta, Sigma = Sigma, 
                                  error_dist = error_dist,
                                  error_parms = error_parms, 
                                  marginal_dists = marginal_dists,
                                  parameters = parameters, 
                                  MissProp = MissProp, 
                                  MissMech = "None")
      
      linear_combs0 <- popMat %*% Eta[-1]
  
      ## Logit fun
      f_alpha <- function(alpha) {
        linear_combs <- alpha + linear_combs0
        prob_miss <- mean(1 / (1 + exp(-linear_combs)))
      }
      
      ## Objective function for optim
      obj_fun <- function(alpha) {
        (f_alpha(alpha) - MissProp)^2
      }
      
      opt_out <- optim(
        par = 0,          # initial guess for alpha
        fn  = obj_fun
      )
      
      Eta[1] <- opt_out$par
    }
    
  
    # Use found intercept to induce missingness into the dataset
    Dataset <- cbind(y = y, X = X)
    linear_combs <- Eta[1] + Dataset %*% as.matrix(Eta[-1])
    miss_vec <- 1 / (1 + exp(-linear_combs))
    
  } else {
    # If intercept is given, no need to find it
    if (length(Eta) != p + 2) {
      
      if (is.null(marginal_means)) {
        stop("Need true marginal means of predictors to use Probit intercept calculations")
        
      } else {
        true_y_mean <- marginal_means %*% Beta[-1]
        true_means <- c(true_y_mean, marginal_means)
      }
      
      Eta <- c(0, Eta)
      # Generate large dataset for covariance matrix
      popMat <- OldDataGeneration(n = 1e6, # 1e8 amount blew up my PC
                                  Beta = Beta, Sigma = Sigma, 
                                  error_dist = error_dist,
                                  error_parms = error_parms, 
                                  marginal_dists = marginal_dists,
                                  parameters = parameters, 
                                  MissProp = MissProp, 
                                  MissMech = "None")
      sigma_snake <- cov(popMat)
      
      ## Probit fun
      f_alpha <- function(alpha) {
        z <- (alpha + t(Eta[-1]) %*% true_means) / 
          sqrt(1 + as.numeric(t(Eta[-1]) %*% sigma_snake %*% Eta[-1]))
        pnorm(z)
      }
      
      ## Objective function for optim
      obj_fun <- function(alpha) {
        (f_alpha(alpha) - MissProp)^2
      }
      
      opt_out <- optim(
        par = 0,          # initial guess for alpha
        fn  = obj_fun
        # method = "SANN"
      )
      
      Eta[1] <- opt_out$par
    }
    
    # Use found intercept to induce missingness into the dataset
    Dataset <- cbind(y = y, X = X)
    miss_vec <- pnorm(Eta[1] + Dataset %*% as.matrix(Eta[-1]))
  }
  
  # Diff missing probs for each row use that for runif cutoff
  Dataset[runif(n) < miss_vec, MissCol] <- NA

  # Names
  colnames(Dataset) <- c("Y", paste("X", 1:(length(Beta) - 1), sep = ""))
  
  return(list(Dataset, Eta))
}

### Complete Case Cleaning
ArrayCleaning <- function(Array){
  # List to store cleaned matrices
  CCData <- vector("list", dim(BootstrappedData)[3])
  
  # Loop through each slice
  for (i in 1:dim(BootstrappedData)[3]) {
    mat <- BootstrappedData[, , i]
    CCData[[i]] <- mat[complete.cases(mat), , drop = FALSE]
  }
  
  return(CCData)
}

### Resampling functions
Bootstrap.Fun <- function(Data, B){
  BootDatasets <- array(data = NA, dim = c(dim(Data), B))
  for(i in 1:B){
    BootDatasets[, , i] <- Data[sample(nrow(Data), nrow(Data), replace = TRUE), ]
  }
  return(BootDatasets)
}

Jackknife.Fun <- function(Data){
  JackknifeDatasets <- array(data = NA, 
                             dim = c(dim(Data)[1] - 1, dim(Data)[2], dim(Data)[1]))
  for(i in 1:dim(Data)[1]){
    JackknifeDatasets[, , i] <- Data[-i, ]
  }
  return(JackknifeDatasets)
}

### Imputation Function
Single.Imp.Fun <- function(Data, FUN, ...) {
  for (j in 1:dim(Data)[3]) {
    for (i in 1:dim(Data)[2]) {
      if (anyNA(Data[, i, j])) {
        # Current column of interest
        dats <- Data[, i, j]
        
        # Which observations are missing?
        Miss <- is.na(dats)
        
        # Calculate imputed value for column based on observed values
        ImpVal <- FUN(dats, na.rm = TRUE, ...)
        
        # Input imputed value for each missing observation in the column
        Data[Miss, i, j] <- ImpVal
      }
    }
  }
  return(Data)
}

Imp.Fun <- function(Data, KNNsize, ...) {
  # Copy data so to be imputed in several ways
  Mean.Det.Imp.Data <- Mean.Stoc.Imp.Data <- 
    Reg.Det.Imp.Data <- Reg.Stoc.Imp.Data <- 
    KNN.Det.Imp.Data <- KNN.Stoc.Imp.Data <- Data
  for (k in 1:dim(Data)[3]) {
    set.seed(l)
    for (j in 1:dim(Data)[2]) {
      # Checking if the jth column of kth layer has any NAs
      if (anyNA(Data[, j, k])) {
        # Current column of interest
        dats <- Data[, j, k]
        
        # Which observations are missing?
        Miss <- is.na(dats)
        
        # Calculate imputed value for column based on observed values and impute missing values
        # Deterministic 
        Mean.Det.Imp.Data[Miss, j, k] <- mean(dats, na.rm = TRUE, ...)
        
        # Stochastic
        Mean.Stoc.Imp.Data[Miss, j, k] <- mean(dats, na.rm = TRUE, ...) + 
          rnorm(sum(Miss), mean = 0, sd = sd(dats, na.rm = TRUE))
        
        for (i in 1:dim(Data)[1]) {
          # If current row has an NA, impute it
          if(is.na(Data[i, j, k])) {
            # Identify missingness pattern in row
            MissPat <- is.na(Data[i, , k])
            
            # Find other rows that have the value we want to predict and the same observed values
            Needed.Var <- c(j, which(!MissPat))
            
            # Without the row of the value we want to predict, with the observed variable
            Pred.Data <- Data[-i, Needed.Var, k]
            
            # Removing NAs for stuff
            Pred.Data <- Pred.Data[complete.cases(Pred.Data), ]
            
            ### KNN Imputation
            KNN.Train <- as.matrix(Pred.Data[, -1])
            KNN.Test <- matrix(Data[i, which(!MissPat), k], nrow = 1)
            
            # Impute NA using KNN
            KNN.Det.Imp.Data[i, j, k] <- knn.reg(train = KNN.Train, 
                                                 test = KNN.Test, 
                                                 y = Pred.Data[, 1], k = KNNsize)[["pred"]]
            
            # Find closest observations to value of interest
            KNN_rslt <- get.knn(rbind(KNN.Test, KNN.Train), k = KNNsize)
            
            # 3 Closest observations and their distances
            CloseObs <- KNN_rslt[["nn.index"]][1, ] - 1
            CloseWeights <- KNN_rslt[["nn.dist"]][1, ]
            
            # Picking one randomly, weights inversely proportional to distance
            Ob <- sample(CloseObs, size = 1, prob = 1/CloseWeights)
            KNN.Stoc.Imp.Data[i, j, k] <- Pred.Data[Ob, 1]
            
            ### Regression Imputation
            NewData <- rbind.data.frame(KNN.Test)
            colnames(NewData) <- paste("X", 1:(ncol(Pred.Data)- 1), sep = "")
            colnames(Pred.Data) <- c("Y", paste("X", 1:(ncol(Pred.Data) - 1), sep = "")) 
            Mod <- lm(Y ~ . , data = as.data.frame(Pred.Data))
            
            # Impute Value
            Reg.Det.Imp.Data[i, j, k] <- predict(Mod, newdata = NewData)
            Reg.Stoc.Imp.Data[i, j, k] <- predict(Mod, newdata = NewData) + 
              rnorm(1, 0, sd = summary(Mod)[["sigma"]])
          }
        }
      }
    }
    if (k %% 5 == 0) {
      cat(paste0(k, "th", collapse = ""), "resampled dataset imputed \n")
    }
  }
  return(list(KNN.Det.Imp.Data = KNN.Det.Imp.Data,
              KNN.Stoc.Imp.Data = KNN.Stoc.Imp.Data,
              Reg.Det.Imp.Data = Reg.Det.Imp.Data,
              Reg.Stoc.Imp.Data = Reg.Stoc.Imp.Data,
              Mean.Det.Imp.Data = Mean.Det.Imp.Data,
              Mean.Stoc.Imp.Data = Mean.Stoc.Imp.Data))
} # There's a seed set at each iter to match parallel

Imp.Parallel.Fun <- function(Data, KNNsize, ...) {
  # Function to process one slice
  SliceProcessor <- function(k) {
    # Print progress (takes a second)
    message("Processing slice:", k)
    
    # Initialize variables for this slice
    Mean.Det.Imp.Data <- Mean.Stoc.Imp.Data <- 
      Reg.Det.Imp.Data <- Reg.Stoc.Imp.Data <- 
      KNN.Det.Imp.Data <- KNN.Stoc.Imp.Data <- Data[, , k]
    
    # Loop through each column to check for NAs
    for (j in 1:dim(Data)[2]) {
      if (anyNA(Data[, j, k])) {
        dats <- Data[, j, k]
        Miss <- is.na(dats)
        
        # Mean imputation
        Mean.Det.Imp.Data[Miss, j] <- mean(dats, na.rm = TRUE)
        Mean.Stoc.Imp.Data[Miss, j] <- mean(dats, na.rm = TRUE) + 
          rnorm(sum(Miss), mean = 0, sd = var(dats, na.rm = TRUE))
        
        # Loop through rows for more detailed imputation
        for (i in 1:dim(Data)[1]) {
          if (is.na(Data[i, j, k])) {
            MissPat <- is.na(Data[i, , k])
            Needed.Var <- c(j, which(!MissPat))
            Pred.Data <- Data[-i, Needed.Var, k]
            Pred.Data <- Pred.Data[complete.cases(Pred.Data), ]
            
            if (nrow(Pred.Data) > 0) {
              # KNN-based imputation
              KNN.Train <- as.matrix(Pred.Data[, -1])
              KNN.Test <- matrix(Data[i, which(!MissPat), k], nrow = 1)
              
              KNN.Det.Imp.Data[i, j] <- knn.reg(train = KNN.Train, test = KNN.Test, 
                                                y = Pred.Data[, 1], k = KNNsize)[["pred"]]
              
              KNN_rslt <- get.knn(rbind(KNN.Test, KNN.Train), k = KNNsize)
              CloseObs <- KNN_rslt[["nn.index"]][1, ] - 1
              CloseWeights <- KNN_rslt[["nn.dist"]][1, ]
              
              # Weighted random sampling for stochastic KNN imputation
              Ob <- sample(CloseObs, size = 1, prob = 1/CloseWeights)
              KNN.Stoc.Imp.Data[i, j] <- Pred.Data[Ob, 1]
              
              # Regression-based imputation
              NewData <- data.frame(KNN.Test)
              colnames(NewData) <- paste("X", 1:(ncol(Pred.Data)-1), sep = "")
              colnames(Pred.Data) <- c("Y", paste("X", 1:(ncol(Pred.Data) - 1), sep = "")) 
              Mod <- lm(Y ~ ., data = as.data.frame(Pred.Data))
              
              Reg.Det.Imp.Data[i, j] <- predict(Mod, newdata = NewData)
              Reg.Stoc.Imp.Data[i, j] <- predict(Mod, newdata = NewData) + 
                rnorm(1, 0, sd = summary(Mod)[["sigma"]])
            }
          }
        }
      }
    }
    
    # Return a list of the 6 imputation slices
    list(
      KNN.Det.Imp.Data = KNN.Det.Imp.Data,
      KNN.Stoc.Imp.Data = KNN.Stoc.Imp.Data,
      Reg.Det.Imp.Data = Reg.Det.Imp.Data,
      Reg.Stoc.Imp.Data = Reg.Stoc.Imp.Data,
      Mean.Det.Imp.Data = Mean.Det.Imp.Data,
      Mean.Stoc.Imp.Data = Mean.Stoc.Imp.Data
    )
  }
  
  # Parallelize the SliceProcessor function
  plan(multisession, workers = parallel::detectCores() - 1)
  
  # Run SliceProcessor across slices using future_lapply
  results <- future_lapply(1:dim(Data)[3], function(k) {
    set.seed(l)  # Set seed for each worker
    SliceProcessor(k)
  })
  
  # Combine results into one list of 6 arrays
  Mean.Det.Imp.Data <- Mean.Stoc.Imp.Data <- 
    Reg.Det.Imp.Data <- Reg.Stoc.Imp.Data <- 
    KNN.Det.Imp.Data <- KNN.Stoc.Imp.Data <- array(data = NA, dim = dim(Data))
  
  for (i in 1:length(results)) {
    KNN.Det.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[1]]]
    KNN.Stoc.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[2]]]
    Reg.Det.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[3]]]
    Reg.Stoc.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[4]]]
    Mean.Det.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[5]]]
    Mean.Stoc.Imp.Data[, , i] <- results[[i]][[names(results[[1]])[6]]]
  }
  
  combined_results <- list(
    KNN.Det.Imp.Data = KNN.Det.Imp.Data,
    KNN.Stoc.Imp.Data = KNN.Stoc.Imp.Data,
    Reg.Det.Imp.Data = Reg.Det.Imp.Data,
    Reg.Stoc.Imp.Data = Reg.Stoc.Imp.Data,
    Mean.Det.Imp.Data = Mean.Det.Imp.Data,
    Mean.Stoc.Imp.Data = Mean.Stoc.Imp.Data
  )
  
  # Return combined results
  return(combined_results)
}
# Both above have an l seed for the LOO method script

# When No Resampling
Imp.NOBS.Fun <- function(Data, KNNsize, ...) {
  # Copy data so to be imputed in several ways
  Mean.Det.Imp.Data <- Mean.Stoc.Imp.Data <- 
    Reg.Det.Imp.Data <- Reg.Stoc.Imp.Data <- 
    KNN.Det.Imp.Data <- KNN.Stoc.Imp.Data <- Data
  for (j in 1:dim(Data)[2]) {
    # Checking if the jth column of kth layer has any NAs
    if (anyNA(Data[, j])) {
      # Current column of interest
      dats <- Data[, j]
      
      # Which observations are missing?
      Miss <- is.na(dats)
      
      # Calculate imputed value for column based on observed values and impute missing values
      # Deterministic 
      Mean.Det.Imp.Data[Miss, j] <- mean(dats, na.rm = TRUE, ...)
      
      # Stochastic
      Mean.Stoc.Imp.Data[Miss, j] <- mean(dats, na.rm = TRUE, ...) + 
        rnorm(sum(Miss), mean = 0, sd = var(dats, na.rm = TRUE))
      
      for (i in 1:dim(Data)[1]) {
        # If current row has an NA, impute it
        if(is.na(Data[i, j])) {
          # Identify missingness pattern in row
          MissPat <- is.na(Data[i, ])
          
          # Find other rows that have the value we want to predict and the same observed values
          Needed.Var <- c(j, which(!MissPat))
          
          # Without the row of the value we want to predict, with the observed variable
          Pred.Data <- Data[-i, Needed.Var]
          
          # Removing NAs for stuff
          Pred.Data <- Pred.Data[complete.cases(Pred.Data), ]
          
          ### KNN Imputation
          KNN.Train <- as.matrix(Pred.Data[, -1])
          KNN.Test <- matrix(Data[i, which(!MissPat)], nrow = 1)
          
          # Impute NA using KNN
          KNN.Det.Imp.Data[i, j] <- knn.reg(train = KNN.Train, test = KNN.Test, 
                                            y = Pred.Data[, 1], k = KNNsize)[["pred"]]
          
          # Find closest observations to value of interest
          KNN_rslt <- get.knn(rbind(KNN.Test, KNN.Train), k = KNNsize)
          
          # 3 Closest observations and their distances
          CloseObs <- KNN_rslt[["nn.index"]][1, ] - 1
          CloseWeights <- KNN_rslt[["nn.dist"]][1, ]
          
          # Picking one randomly, weights inversely proportional to distance
          Ob <- sample(CloseObs, size = 1, prob = 1/CloseWeights)
          KNN.Stoc.Imp.Data[i, j] <- Pred.Data[Ob, 1]
          
          ### Regression Imputation
          NewData <- rbind.data.frame(KNN.Test)
          colnames(NewData) <- paste("X", 1:(ncol(Pred.Data)- 1), sep = "")
          colnames(Pred.Data) <- c("Y", paste("X", 1:(ncol(Pred.Data) - 1), sep = "")) 
          Mod <- lm(Y ~ . , data = as.data.frame(Pred.Data))
          
          # Impute Value
          Reg.Det.Imp.Data[i, j] <- predict(Mod, newdata = NewData)
          Reg.Stoc.Imp.Data[i, j] <- predict(Mod, newdata = NewData) + 
            rnorm(1, 0, sd = summary(Mod)[["sigma"]])
        }
      }
    }
  }
  return(list(KNN.Det.Imp.Data = KNN.Det.Imp.Data,
              KNN.Stoc.Imp.Data = KNN.Stoc.Imp.Data,
              Reg.Det.Imp.Data = Reg.Det.Imp.Data,
              Reg.Stoc.Imp.Data = Reg.Stoc.Imp.Data,
              Mean.Det.Imp.Data = Mean.Det.Imp.Data,
              Mean.Stoc.Imp.Data = Mean.Stoc.Imp.Data))
}

### Model Selection Functions
MS.Fun <- function(Data, Beta, criteria = "AIC"){
  # Number of observations in each bootstrap dataset, assuming each element slice is the same height
  nsamp <- dim(Data[[1]])[1]
  
  # Number of parameters (+ intercept) assuming that each list element and slice is the same width
  BetaL <- dim(Data[[1]])[2]
  
  # Number of bootstrap datasets assuming each list element has the same amount of slices
  BootSamps <- dim(Data[[1]])[3]
  
  # Initialize progress bar
  prog <- progressor(steps = length(Data)*BootSamps)
  
  # Create object for selected coefficients for each imputation method and each bootstrap dataset
  Ridge.Coef.Array <- LASSO.Coef.Array <- array(NA, dim = c(BootSamps, (BetaL - 1), length(Data)))
  BestSub.Reg.Sel.List <- Forward.Reg.Sel.List <- lapply(1:length(Data), \(x) list())
  
  # Create objects to store coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef.Array <- SCoef.Array <- CCoef.Array <- FCoef.Array <- array(NA, dim = c(BootSamps, BetaL, length(Data)))
  BS.CoefSE.Array <- SCoefSE.Array <- CCoefSE.Array <- FCoefSE.Array <- array(NA, dim = c(BootSamps, BetaL, length(Data)))
  BS.T.Val.Array <- ST.Val.Array <- CT.Val.Array <- FT.Val.Array <- array(NA, dim = c(BootSamps, BetaL, length(Data)))
  BS.Res.SE <- S.Res.SE <- C.Res.SE <- F.Res.SE <- matrix(NA, nrow = BootSamps, ncol = length(Data))
  
  for(j in 1:length(Data)){
    for(i in 1:BootSamps){
      # Copy in objects for ith array
      y.vec <- Data[[j]][, 1, i]
      x.mat <- Data[[j]][, -1, i]
      
      ### Penalized Regression Methods
      ## LASSO
      # Cross validation to calculate lambda
      LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                            standardize = TRUE, thresh = 1e-10)
      # Fitting Model
      LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                          standardize = TRUE, thresh = 1e-10)
      # Pulling coefficients
      LASSO.Coef.Array[i, , j] <- as.matrix(LASSO.mod$beta)
      
      ## Ridge
      # Cross validation to calculate lambda
      Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
      # Fitting Model
      Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                          standardize = TRUE, thresh = 1e-10)
      # Pulling coefficients
      Ridge.Coef.Array[i, , j] <- as.matrix(Ridge.mod$beta)
      
      ### Sequential Methods
      da.data <- as.data.frame(Data[[j]][, , i])
      colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
      null.model <- lm(Y ~ 1, data = da.data)
      
      ## Forward Selection
      Mod.Select.For <- step(object = null.model, 
                             scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                             direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
      
      # Noting which regressors were chosen in the selected model
      Forward.Reg.Sel.List[[j]][[i]] <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
      
      # Pulling summary from selected model
      lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
      
      # Checking which index the variables chosen correspond to
      ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
      
      # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
      SCoef.Array[i, c(1, ind + 1), j] <- lm.select$coefficients[, 1]
      SCoefSE.Array[i, c(1, ind + 1), j] <- lm.select$coefficients[, 2]
      ST.Val.Array[i, c(1, ind + 1), j] <- lm.select$coefficients[, 3]
      S.Res.SE[i, j] <- lm.select$sigma
      
      # Pulling summary from"correct" model
      # Checking for which regressors should be in the model
      Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
      lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
      ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
      
      # Saving relevant statistics
      CCoef.Array[i, c(1, ind + 1), j] <- lm.correct$coefficients[, 1]
      CCoefSE.Array[i, c(1, ind + 1), j] <- lm.correct$coefficients[, 2]
      CT.Val.Array[i, c(1, ind + 1), j] <- lm.correct$coefficients[, 3]
      C.Res.SE[i, j] <- lm.correct$sigma
      
      # Pulling summary from full model
      lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
      ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
      
      # Saving relevant statistics
      FCoef.Array[i, c(1, ind + 1), j] <- lm.full$coefficients[, 1]
      FCoefSE.Array[i, c(1, ind + 1), j] <- lm.full$coefficients[, 2]
      FT.Val.Array[i, c(1, ind + 1), j] <- lm.full$coefficients[, 3]
      F.Res.SE[i, j] <- lm.full$sigma
      
      ## Backward Selection
      # Add later if wanted
      
      ### Best Subsets
      BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
      ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
      
      if(criteria == "AIC"){
        # Adding AIC, using formula in step() function
        n <- nrow(da.data)
        ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
        
        # Which regressors where chosen?
        BestSub.Reg.Sel.List[[j]][[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
        
      } else if (criteria == "BIC") {
        # Which regressors where chosen?
        BestSub.Reg.Sel.List[[j]][[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
        
      } else {
        print("Criteria not supported")
      }
      
      # Pulling summary from selected model
      Vars <- strsplit(BestSub.Reg.Sel.List[[j]][[i]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
      LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
      
      # Not indexing here, hoping that won't make me cry later - should be in correct order...
      # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
      BS.Coef.Array[i, , j] <- LM.Select$coefficients[, 1]
      BS.CoefSE.Array[i, , j] <- LM.Select$coefficients[, 2]
      BS.T.Val.Array[i, , j] <- LM.Select$coefficients[, 3]
      BS.Res.SE[i, j] <- LM.Select$sigma
      
      # Update progress bar
      prog()
    }
  }
  
  # Output naming
  ImpNames <- names(Data)
  dimnames(LASSO.Coef.Array) <- dimnames(Ridge.Coef.Array) <- 
    list(NULL, paste("Beta", 1:(BetaL - 1), sep = ""), ImpNames)
  dimnames(BS.Coef.Array) <- dimnames(SCoef.Array) <- dimnames(CCoef.Array) <- 
    dimnames(FCoef.Array) <- list(NULL, c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.CoefSE.Array) <- dimnames(SCoefSE.Array) <- dimnames(CCoefSE.Array) <- dimnames(FCoefSE.Array) <- 
    list(NULL, c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.T.Val.Array) <- dimnames(ST.Val.Array) <- dimnames(CT.Val.Array) <- dimnames(FT.Val.Array) <- 
    list(NULL, c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  names(Forward.Reg.Sel.List) <- names(BestSub.Reg.Sel.List) <- ImpNames
  for (method in ImpNames) {
    names(Forward.Reg.Sel.List[[method]]) <- names(BestSub.Reg.Sel.List[[method]]) <- paste0("BootData", 1:BootSamps)
  }
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef.Array, Ridge = Ridge.Coef.Array)
  ForwardSelection <- list(Coefficients = SCoef.Array, StandardErrors = SCoefSE.Array, TStats = ST.Val.Array, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef.Array, StandardErrors = BS.CoefSE.Array, TStats = BS.T.Val.Array, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef.Array, StandardErrors = CCoefSE.Array, TStats = CT.Val.Array, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef.Array, StandardErrors = FCoefSE.Array, TStats = FT.Val.Array, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel.List, ForwardSelection = Forward.Reg.Sel.List)
  
  Output <- list(PenalizedRegression = PenalizedRegression,
                 ForwardSelection = ForwardSelection,
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

MS.Stack.Fun <- function(Data, Beta, criteria = "AIC"){
  # Number of observations in each bootstrap dataset, assuming each element slice is the same height
  nsamp <- dim(Data[[1]])[1]
  
  # Number of parameters (+ intercept) assuming that each list element and slice is the same width
  BetaL <- dim(Data[[1]])[2]
  
  # Initialize progress bar
  prog <- progressor(steps = length(Data))
  
  # Create object for selected coefficients for each imputation method and each bootstrap dataset
  Ridge.Coef.Matrix <- LASSO.Coef.Matrix <- matrix(NA, nrow = (BetaL - 1), ncol = length(Data))
  BestSub.Reg.Sel.List <- Forward.Reg.Sel.List <- vector("list", length(Data))
  
  # Create objects to store coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef.Matrix <- SCoef.Matrix <- CCoef.Matrix <- FCoef.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.CoefSE.Matrix <- SCoefSE.Matrix <- CCoefSE.Matrix <- FCoefSE.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.T.Val.Matrix <- ST.Val.Matrix <- CT.Val.Matrix <- FT.Val.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.Res.SE <- S.Res.SE <- C.Res.SE <- F.Res.SE <- vector("double", length(Data))
  
  for(i in 1:length(Data)){
    # Stack all resampled imputed datasets for each imputation method
    BigMat <- do.call(rbind, lapply(1:dim(Data[[i]])[3], function(k) Data[[i]][, , k]))
    
    # Copy in X and Y objects
    y.vec <- BigMat[, 1]
    x.mat <- BigMat[, -1]
    
    ### Penalized Regression Methods
    ## LASSO
    # Cross validation to calculate lambda
    LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                          standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    LASSO.Coef.Matrix[, i] <- as.matrix(LASSO.mod$beta)
    
    ## Ridge
    # Cross validation to calculate lambda
    Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    Ridge.Coef.Matrix[, i] <- as.matrix(Ridge.mod$beta)
    
    ### Sequential Methods
    da.data <- as.data.frame(BigMat)
    colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
    null.model <- lm(Y ~ 1, data = da.data)
    
    ## Forward Selection
    Mod.Select.For <- step(object = null.model, 
                           scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                           direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
    
    # Noting which regressors were chosen in the selected model
    Forward.Reg.Sel.List[[i]] <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
    
    # Pulling summary from selected model
    lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
    
    # Checking which index the variables chosen correspond to
    ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    SCoef.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 1]
    SCoefSE.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 2]
    ST.Val.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 3]
    S.Res.SE[i] <- lm.select$sigma
    
    # Pulling summary from"correct" model
    # Checking for which regressors should be in the model
    Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
    lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
    ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    CCoef.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 1]
    CCoefSE.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 2]
    CT.Val.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 3]
    C.Res.SE[i] <- lm.correct$sigma
    
    # Pulling summary from full model
    lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
    ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    FCoef.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 1]
    FCoefSE.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 2]
    FT.Val.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 3]
    F.Res.SE[i] <- lm.full$sigma
    
    ## Backward Selection
    # Add later if wanted
    
    ### Best Subsets
    BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
    ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
    
    if(criteria == "AIC"){
      # Adding AIC, using formula in step() function
      n <- nrow(da.data)
      ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
      
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
      
    } else if (criteria == "BIC") {
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
      
    } else {
      print("Criteria not supported")
    }
    
    # Pulling summary from selected model
    Vars <- strsplit(BestSub.Reg.Sel.List[[i]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
    LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
    
    # Not indexing here, hoping that won't make me cry later - should be in correct order...
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    BS.Coef.Matrix[, i] <- LM.Select$coefficients[, 1]
    BS.CoefSE.Matrix[, i] <- LM.Select$coefficients[, 2]
    BS.T.Val.Matrix[, i] <- LM.Select$coefficients[, 3]
    BS.Res.SE[i] <- LM.Select$sigma
    
    # Update progress bar
    prog()
  }
  
  # Output naming
  ImpNames <- names(Data)
  dimnames(LASSO.Coef.Matrix) <- dimnames(Ridge.Coef.Matrix) <- 
    list(paste("Beta", 1:(BetaL - 1), sep = ""), ImpNames)
  dimnames(BS.Coef.Matrix) <- dimnames(SCoef.Matrix) <- dimnames(CCoef.Matrix) <- 
    dimnames(FCoef.Matrix) <- list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.CoefSE.Matrix) <- dimnames(SCoefSE.Matrix) <- dimnames(CCoefSE.Matrix) <- dimnames(FCoefSE.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.T.Val.Matrix) <- dimnames(ST.Val.Matrix) <- dimnames(CT.Val.Matrix) <- dimnames(FT.Val.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  names(Forward.Reg.Sel.List) <- names(BestSub.Reg.Sel.List) <- ImpNames
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef.Matrix, Ridge = Ridge.Coef.Matrix)
  ForwardSelection <- list(Coefficients = SCoef.Matrix, StandardErrors = SCoefSE.Matrix, TStats = ST.Val.Matrix, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef.Matrix, StandardErrors = BS.CoefSE.Matrix, TStats = BS.T.Val.Matrix, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef.Matrix, StandardErrors = CCoefSE.Matrix, TStats = CT.Val.Matrix, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef.Matrix, StandardErrors = FCoefSE.Matrix, TStats = FT.Val.Matrix, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel.List, ForwardSelection = Forward.Reg.Sel.List)
  
  Output <- list(PenalizedRegression = PenalizedRegression, 
                 ForwardSelection = ForwardSelection, 
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

# Handling Complete Cases
MS.CC.Fun <- function(Data, Beta, criteria = "AIC"){
  # Number of parameters (+ intercept) assuming that each list element matrix has the same width
  BetaL <- dim(Data[[1]])[2]
  
  # Number of bootstrap datasets
  BootSamps <- length(Data)
  
  # Initialize progress bar
  prog <- progressor(steps = BootSamps)
  
  # Create object for selected coefficients for each imputation method and each bootstrap dataset
  Ridge.Coef.Matrix <- LASSO.Coef.Matrix <- matrix(NA, nrow = (BetaL - 1), ncol = length(Data))
  BestSub.Reg.Sel.List <- Forward.Reg.Sel.List <- vector("list", length(Data))
  
  # Create objects to store coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef.Matrix <- SCoef.Matrix <- CCoef.Matrix <- FCoef.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.CoefSE.Matrix <- SCoefSE.Matrix <- CCoefSE.Matrix <- FCoefSE.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.T.Val.Matrix <- ST.Val.Matrix <- CT.Val.Matrix <- FT.Val.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.Res.SE <- S.Res.SE <- C.Res.SE <- F.Res.SE <- vector("double", length(Data))
  
  for(i in 1:BootSamps){
    # Copy in objects for ith list element
    y.vec <- Data[[i]][, 1]
    x.mat <- Data[[i]][, -1]
    
    ### Penalized Regression Methods
    ## LASSO
    # Cross validation to calculate lambda
    LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                          standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    LASSO.Coef.Matrix[, i] <- as.matrix(LASSO.mod$beta)
    
    ## Ridge
    # Cross validation to calculate lambda
    Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    Ridge.Coef.Matrix[, i] <- as.matrix(Ridge.mod$beta)
    
    ### Sequential Methods
    da.data <- as.data.frame(Data[[i]])
    colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
    null.model <- lm(Y ~ 1, data = da.data)
    
    ## Forward Selection
    Mod.Select.For <- step(object = null.model, 
                           scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                           direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
    
    # Noting which regressors were chosen in the selected model
    Forward.Reg.Sel.List[[i]] <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
    
    # Pulling summary from selected model
    lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
    
    # Checking which index the variables chosen correspond to
    ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    SCoef.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 1]
    SCoefSE.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 2]
    ST.Val.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 3]
    S.Res.SE[i] <- lm.select$sigma
    
    # Pulling summary from"correct" model
    # Checking for which regressors should be in the model
    Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
    lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
    ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    CCoef.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 1]
    CCoefSE.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 2]
    CT.Val.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 3]
    C.Res.SE[i] <- lm.correct$sigma
    
    # Pulling summary from full model
    lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
    ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    FCoef.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 1]
    FCoefSE.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 2]
    FT.Val.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 3]
    F.Res.SE[i] <- lm.full$sigma
    
    ## Backward Selection
    # Add later if wanted
    
    ### Best Subsets
    BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
    ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
    
    if(criteria == "AIC"){
      # Adding AIC, using formula in step() function
      n <- nrow(da.data)
      ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
      
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
      
    } else if (criteria == "BIC") {
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
      
    } else {
      print("Criteria not supported")
    }
    
    # Pulling summary from selected model
    Vars <- strsplit(BestSub.Reg.Sel.List[[i]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
    LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
    
    # Not indexing here, hoping that won't make me cry later - should be in correct order...
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    BS.Coef.Matrix[, i] <- LM.Select$coefficients[, 1]
    BS.CoefSE.Matrix[, i] <- LM.Select$coefficients[, 2]
    BS.T.Val.Matrix[, i] <- LM.Select$coefficients[, 3]
    BS.Res.SE[i] <- LM.Select$sigma
    
    # Update progress bar
    prog()
  }
  
  # Output naming
  dimnames(LASSO.Coef.Matrix) <- dimnames(Ridge.Coef.Matrix) <- 
    list(paste("Beta", 1:(BetaL - 1), sep = ""), paste("DS", 1:length(Data), sep = ""))
  dimnames(BS.Coef.Matrix) <- dimnames(SCoef.Matrix) <- dimnames(CCoef.Matrix) <- 
    dimnames(FCoef.Matrix) <- list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), paste("DS", 1:length(Data), sep = ""))
  dimnames(BS.CoefSE.Matrix) <- dimnames(SCoefSE.Matrix) <- dimnames(CCoefSE.Matrix) <- dimnames(FCoefSE.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), paste("DS", 1:length(Data), sep = ""))
  dimnames(BS.T.Val.Matrix) <- dimnames(ST.Val.Matrix) <- dimnames(CT.Val.Matrix) <- dimnames(FT.Val.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), paste("DS", 1:length(Data), sep = ""))
  names(Forward.Reg.Sel.List) <- names(BestSub.Reg.Sel.List) <- paste("DS", 1:length(Data), sep = "")
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef.Matrix, Ridge = Ridge.Coef.Matrix)
  ForwardSelection <- list(Coefficients = SCoef.Matrix, StandardErrors = SCoefSE.Matrix, TStats = ST.Val.Matrix, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef.Matrix, StandardErrors = BS.CoefSE.Matrix, TStats = BS.T.Val.Matrix, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef.Matrix, StandardErrors = CCoefSE.Matrix, TStats = CT.Val.Matrix, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef.Matrix, StandardErrors = FCoefSE.Matrix, TStats = FT.Val.Matrix, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel.List, ForwardSelection = Forward.Reg.Sel.List)
  
  Output <- list(PenalizedRegression = PenalizedRegression, 
                 ForwardSelection = ForwardSelection, 
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

MS.CC.Stack.Fun <- function(Data, Beta, criteria = "AIC"){
  # Stack all resampled imputed datasets for each imputation method
  BigMat <- do.call(rbind, Data)
  
  # Copy in X and Y objects
  y.vec <- BigMat[, 1]
  x.mat <- BigMat[, -1]
  
  ### Penalized Regression Methods
  ## LASSO
  # Cross validation to calculate lambda
  LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                        standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  LASSO.Coef <- as.matrix(LASSO.mod$beta)
  
  ## Ridge
  # Cross validation to calculate lambda
  Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  Ridge.Coef <- as.matrix(Ridge.mod$beta)
  
  ### Sequential Methods
  da.data <- as.data.frame(BigMat)
  colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
  null.model <- lm(Y ~ 1, data = da.data)
  
  ## Forward Selection
  Mod.Select.For <- step(object = null.model, 
                         scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                         direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
  
  # Noting which regressors were chosen in the selected model
  Forward.Reg.Sel <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
  
  # Pulling summary from selected model
  lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
  
  # Checking which index the variables chosen correspond to
  ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
  SCoef <- lm.select$coefficients[, 1]
  SCoefSE <- lm.select$coefficients[, 2]
  ST.Val <- lm.select$coefficients[, 3]
  S.Res.SE <- lm.select$sigma
  
  # Pulling summary from"correct" model
  # Checking for which regressors should be in the model
  Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
  lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
  ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving relevant statistics
  CCoef <- lm.correct$coefficients[, 1]
  CCoefSE <- lm.correct$coefficients[, 2]
  CT.Val <- lm.correct$coefficients[, 3]
  C.Res.SE <- lm.correct$sigma
  
  # Pulling summary from full model
  lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
  ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving relevant statistics
  FCoef <- lm.full$coefficients[, 1]
  FCoefSE <- lm.full$coefficients[, 2]
  FT.Val <- lm.full$coefficients[, 3]
  F.Res.SE <- lm.full$sigma
  
  ## Backward Selection
  # Add later if wanted
  
  ### Best Subsets
  BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
  ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
  
  if(criteria == "AIC"){
    # Adding AIC, using formula in step() function
    n <- nrow(da.data)
    ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
    
    # Which regressors where chosen?
    BestSub.Reg.Sel <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
    
  } else if (criteria == "BIC") {
    # Which regressors where chosen?
    BestSub.Reg.Sel <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
    
  } else {
    print("Criteria not supported")
  }
  
  # Pulling summary from selected model
  Vars <- strsplit(BestSub.Reg.Sel, "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
  LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
  
  # Not indexing here, hoping that won't make me cry later - should be in correct order...
  # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef <- LM.Select$coefficients[, 1]
  BS.CoefSE <- LM.Select$coefficients[, 2]
  BS.T.Val <- LM.Select$coefficients[, 3]
  BS.Res.SE <- LM.Select$sigma
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef, Ridge = Ridge.Coef)
  ForwardSelection <- list(Coefficients = SCoef, StandardErrors = SCoefSE, TStats = ST.Val, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef, StandardErrors = BS.CoefSE, TStats = BS.T.Val, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef, StandardErrors = CCoefSE, TStats = CT.Val, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef, StandardErrors = FCoefSE, TStats = FT.Val, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel, ForwardSelection = Forward.Reg.Sel)
  
  Output <- list(PenalizedRegression = PenalizedRegression, 
                 ForwardSelection = ForwardSelection, 
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

# When No Resampling
MS.NOBS.Fun <- function(Data, Beta, criteria = "AIC"){
  # Number of observations in each bootstrap dataset, assuming each element is the same height
  nsamp <- dim(Data[[1]])[1]
  
  # Number of parameters (+ intercept) assuming that each list element is the same width
  BetaL <- dim(Data[[1]])[2]
  
  # Initialize progress bar
  prog <- progressor(steps = length(Data))
  
  # Create object for selected coefficients for each imputation method and each bootstrap dataset
  Ridge.Coef.Matrix <- LASSO.Coef.Matrix <- matrix(NA, nrow = (BetaL - 1), ncol = length(Data))
  BestSub.Reg.Sel.List <- Forward.Reg.Sel.List <- vector("list", length(Data))
  
  # Create objects to store coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef.Matrix <- SCoef.Matrix <- CCoef.Matrix <- FCoef.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.CoefSE.Matrix <- SCoefSE.Matrix <- CCoefSE.Matrix <- FCoefSE.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.T.Val.Matrix <- ST.Val.Matrix <- CT.Val.Matrix <- FT.Val.Matrix <- matrix(NA, nrow = BetaL, ncol = length(Data))
  BS.Res.SE <- S.Res.SE <- C.Res.SE <- F.Res.SE <- vector("double", length(Data))
  
  for(i in 1:length(Data)){
    # Copy in X and Y objects
    y.vec <- Data[[i]][, 1]
    x.mat <- Data[[i]][, -1]
    
    ### Penalized Regression Methods
    ## LASSO
    # Cross validation to calculate lambda
    LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                          standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    LASSO.Coef.Matrix[, i] <- as.matrix(LASSO.mod$beta)
    
    ## Ridge
    # Cross validation to calculate lambda
    Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    # Fitting Model
    Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    Ridge.Coef.Matrix[, i] <- as.matrix(Ridge.mod$beta)
    
    ### Sequential Methods
    da.data <- as.data.frame(Data[[i]])
    colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
    null.model <- lm(Y ~ 1, data = da.data)
    
    ## Forward Selection
    Mod.Select.For <- step(object = null.model, 
                           scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                           direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
    
    # Noting which regressors were chosen in the selected model
    Forward.Reg.Sel.List[[i]] <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
    
    # Pulling summary from selected model
    lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
    
    # Checking which index the variables chosen correspond to
    ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    SCoef.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 1]
    SCoefSE.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 2]
    ST.Val.Matrix[c(1, ind + 1), i] <- lm.select$coefficients[, 3]
    S.Res.SE[i] <- lm.select$sigma
    
    # Pulling summary from"correct" model
    # Checking for which regressors should be in the model
    Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
    lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
    ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    CCoef.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 1]
    CCoefSE.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 2]
    CT.Val.Matrix[c(1, ind + 1), i] <- lm.correct$coefficients[, 3]
    C.Res.SE[i] <- lm.correct$sigma
    
    # Pulling summary from full model
    lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
    ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
    
    # Saving relevant statistics
    FCoef.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 1]
    FCoefSE.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 2]
    FT.Val.Matrix[c(1, ind + 1), i] <- lm.full$coefficients[, 3]
    F.Res.SE[i] <- lm.full$sigma
    
    ## Backward Selection
    # Add later if wanted
    
    ### Best Subsets
    BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
    ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
    
    if(criteria == "AIC"){
      # Adding AIC, using formula in step() function
      n <- nrow(da.data)
      ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
      
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
      
    } else if (criteria == "BIC") {
      # Which regressors where chosen?
      BestSub.Reg.Sel.List[[i]] <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
      
    } else {
      print("Criteria not supported")
    }
    
    # Pulling summary from selected model
    Vars <- strsplit(BestSub.Reg.Sel.List[[i]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
    LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
    
    # Not indexing here, hoping that won't make me cry later - should be in correct order...
    # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
    BS.Coef.Matrix[, i] <- LM.Select$coefficients[, 1]
    BS.CoefSE.Matrix[, i] <- LM.Select$coefficients[, 2]
    BS.T.Val.Matrix[, i] <- LM.Select$coefficients[, 3]
    BS.Res.SE[i] <- LM.Select$sigma
    
    # Update progress bar
    prog()
  }
  
  # Output naming
  ImpNames <- names(Data)
  dimnames(LASSO.Coef.Matrix) <- dimnames(Ridge.Coef.Matrix) <- 
    list(paste("Beta", 1:(BetaL - 1), sep = ""), ImpNames)
  dimnames(BS.Coef.Matrix) <- dimnames(SCoef.Matrix) <- dimnames(CCoef.Matrix) <- 
    dimnames(FCoef.Matrix) <- list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.CoefSE.Matrix) <- dimnames(SCoefSE.Matrix) <- dimnames(CCoefSE.Matrix) <- dimnames(FCoefSE.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  dimnames(BS.T.Val.Matrix) <- dimnames(ST.Val.Matrix) <- dimnames(CT.Val.Matrix) <- dimnames(FT.Val.Matrix) <- 
    list(c("Intercept", paste("Beta", 1:(BetaL - 1), sep = "")), ImpNames)
  names(Forward.Reg.Sel.List) <- names(BestSub.Reg.Sel.List) <- ImpNames
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef.Matrix, Ridge = Ridge.Coef.Matrix)
  ForwardSelection <- list(Coefficients = SCoef.Matrix, StandardErrors = SCoefSE.Matrix, TStats = ST.Val.Matrix, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef.Matrix, StandardErrors = BS.CoefSE.Matrix, TStats = BS.T.Val.Matrix, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef.Matrix, StandardErrors = CCoefSE.Matrix, TStats = CT.Val.Matrix, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef.Matrix, StandardErrors = FCoefSE.Matrix, TStats = FT.Val.Matrix, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel.List, ForwardSelection = Forward.Reg.Sel.List)
  
  Output <- list(PenalizedRegression = PenalizedRegression, 
                 ForwardSelection = ForwardSelection, 
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

# When No Resampling or Imputation
Std.MS.Fun <- function(Data, Beta, criteria = "AIC"){
  # Copy in X and Y objects
  y.vec <- Data[, 1]
  x.mat <- Data[, -1]
  
  ### Penalized Regression Methods
  ## LASSO
  # Cross validation to calculate lambda
  LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, 
                        standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  LASSO.Coef <- as.matrix(LASSO.mod$beta)
  
  ## Ridge
  # Cross validation to calculate lambda
  Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  Ridge.Coef <- as.matrix(Ridge.mod$beta)
  
  ### Sequential Methods
  da.data <- as.data.frame(Data)
  colnames(da.data) <- c("Y", paste("X", 1:(ncol(da.data) - 1), sep = ""))
  null.model <- lm(Y ~ 1, data = da.data)
  
  ## Forward Selection
  Mod.Select.For <- step(object = null.model, 
                         scope = as.formula(paste("~", paste(colnames(da.data)[-1], collapse = " + "))), 
                         direction = "forward", trace = 0, k = ifelse(criteria == "AIC", 2, log(nsamp)))
  
  # Noting which regressors were chosen in the selected model
  Forward.Reg.Sel <- paste(sort(attr(Mod.Select.For$terms, "term.labels")), collapse = "")
  
  # Pulling summary from selected model
  lm.select <- summary(lm(formula = Mod.Select.For$call$formula, data = da.data))
  
  # Checking which index the variables chosen correspond to
  ind <- match(attr(Mod.Select.For$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
  SCoef <- lm.select$coefficients[, 1]
  SCoefSE <- lm.select$coefficients[, 2]
  ST.Val <- lm.select$coefficients[, 3]
  S.Res.SE <- lm.select$sigma
  
  # Pulling summary from"correct" model
  # Checking for which regressors should be in the model
  Correct.Regs <- colnames(da.data)[-1][which(Beta[-1] != 0)]
  lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
  ind <- match(attr(lm.correct$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving relevant statistics
  CCoef <- lm.correct$coefficients[, 1]
  CCoefSE <- lm.correct$coefficients[, 2]
  CT.Val <- lm.correct$coefficients[, 3]
  C.Res.SE <- lm.correct$sigma
  
  # Pulling summary from full model
  lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(da.data)[-1], collapse = " + "))), data = da.data))
  ind <- match(attr(lm.full$terms, "term.labels"), colnames(da.data)[-1])
  
  # Saving relevant statistics
  FCoef <- lm.full$coefficients[, 1]
  FCoefSE <- lm.full$coefficients[, 2]
  FT.Val <- lm.full$coefficients[, 3]
  F.Res.SE <- lm.full$sigma
  
  ## Backward Selection
  # Add later if wanted
  
  ### Best Subsets
  BestSubMods <- summary(regsubsets(Y ~ ., data = da.data, method = "exhaustive", nbest = 1, intercept = TRUE))
  ModelSummary <- with(BestSubMods, data.frame(p = rowSums(which), RSS = rss, AdjR2 = adjr2, Cp = cp, BIC = bic))
  
  if(criteria == "AIC"){
    # Adding AIC, using formula in step() function
    n <- nrow(da.data)
    ModelSummary <- with(ModelSummary, data.frame(ModelSummary, AIC = n * log(RSS / n) + 2 * p))
    
    # Which regressors where chosen?
    BestSub.Reg.Sel <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["AIC"]]), ])[-1], collapse = "")
    
  } else if (criteria == "BIC") {
    # Which regressors where chosen?
    BestSub.Reg.Sel <- paste(names(BestSubMods[["which"]][which.min(ModelSummary[["BIC"]]), ])[-1], collapse = "")
    
  } else {
    print("Criteria not supported")
  }
  
  # Pulling summary from selected model
  Vars <- strsplit(BestSub.Reg.Sel, "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
  LM.Select <- summary(lm(formula = as.formula(paste("Y ~", paste(Vars, collapse = "+"))), data = da.data))
  
  # Not indexing here, hoping that won't make me cry later - should be in correct order...
  # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
  BS.Coef <- LM.Select$coefficients[, 1]
  BS.CoefSE <- LM.Select$coefficients[, 2]
  BS.T.Val <- LM.Select$coefficients[, 3]
  BS.Res.SE <- LM.Select$sigma
  
  # Output structuring
  PenalizedRegression <- list(LASS0 = LASSO.Coef, Ridge = Ridge.Coef)
  ForwardSelection <- list(Coefficients = SCoef, StandardErrors = SCoefSE, TStats = ST.Val, ResidualSE = S.Res.SE)
  BestSubSelection <- list(Coefficients = BS.Coef, StandardErrors = BS.CoefSE, TStats = BS.T.Val, ResidualSE = BS.Res.SE)
  CorrectModel <- list(Coefficients = CCoef, StandardErrors = CCoefSE, TStats = CT.Val, ResidualSE = C.Res.SE)
  FullModel <- list(Coefficients = FCoef, StandardErrors = FCoefSE, TStats = FT.Val, ResidualSE = F.Res.SE)
  ChosenModels <- list(BestSubset = BestSub.Reg.Sel, ForwardSelection = Forward.Reg.Sel)
  
  Output <- list(PenalizedRegression = PenalizedRegression, 
                 ForwardSelection = ForwardSelection, 
                 BestSubSelection = BestSubSelection, 
                 CorrectModel = CorrectModel, 
                 FullModel = FullModel, 
                 ChosenModels = ChosenModels)
  
  return(Output)
}

## Stability Selection Functions
Stability.Selection.Fun <- function(Data, pip){
  # Check is proportion is in 0 to 1 interval
  if(pip >= 0 & pip <= 1){
    
    # Assuming each selection method uses the same imputation methods
    ImpMethods <- names(Data$ChosenModels$ForwardSelection)
    
    # Assuming the methods obtained the same amount of results as bootstrap samples
    Bootsamps <- nrow(Data$PenalizedRegression$LASS0)
    
    # Parsing the regressors that were selected in all the imputed datasets
    ModelSelectionMethods <- c("ForwardSelection", "BestSubset")
    Reg.Together.Tables <- Reg.Sep <- list()
    for(j in ModelSelectionMethods){
      for (i in ImpMethods) {
        Reg.Together.Tables[[j]][[i]] <- table(unlist(Data$ChosenModels[[j]][[i]]))
        for(k in 1:Bootsamps){
          Reg.Sep[[j]][[i]][[k]] <- strsplit(Data$ChosenModels[[j]][[i]][[k]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
        }
      }
    }
    
    # Table of counts for each regressor for each imputation method
    Selection.Regs <- Reg.Sep.Tables <- list()
    for(j in ModelSelectionMethods){
      for (i in ImpMethods) {
        Reg.Sep.Tables[[j]][[i]] <- table(unlist(Reg.Sep[[j]][[i]]))
        Selection.Regs[[j]][[i]] <- names(Reg.Sep.Tables[[j]][[i]][Reg.Sep.Tables[[j]][[i]] >= pip * Bootsamps])
      }
    }
  } else {
    stop("Proportion must be between 0 and 1")
  }
  
  return(list(Reg.Sep.Tables = Reg.Sep.Tables, 
              Selection.Regs = Selection.Regs))
}

# Handling Complete Cases
Stability.Selection.CC.Fun <- function(Data, pip){
  # Check is proportion is in 0 to 1 interval
  if(pip >= 0 & pip <= 1){
    
    # Assuming the number of boostrap samples is the same for each MS method
    Bootsamps <- length(Data$ChosenModels$ForwardSelection)
    
    # Parsing the regressors that were selected
    ModelSelectionMethods <- c("ForwardSelection", "BestSubset")
    Selection.Regs <- Reg.Sep.Tables <- Reg.Sep <- list()
    for(j in ModelSelectionMethods){
      for(i in 1:Bootsamps){
        Reg.Sep[[j]][[i]] <- strsplit(Data$ChosenModels[[j]][[i]], "(?<=\\d)(?=[A-Za-z])", perl = TRUE)[[1]]
      }
      Reg.Sep.Tables[[j]] <- table(unlist(Reg.Sep[[j]]))
      Selection.Regs[[j]] <- names(Reg.Sep.Tables[[j]][Reg.Sep.Tables[[j]] >= pip * Bootsamps])
    }
  } else {
    stop("Proportion must be between 0 and 1")
  }
  
  return(list(Reg.Sep.Tables = Reg.Sep.Tables, 
              Selection.Regs = Selection.Regs))
}

# For fitting the models that pop out of the stability selection functions
FitStabilitySelectedModel <- function(Response, Regressors, data) {
  # Fit the linear model
  model <- lm(as.formula(paste(Response, "~", paste(Regressors, collapse = " + "))), 
              data = as.data.frame(data))
  
  # Extract summary statistics
  model_summary <- summary(model)
  
  # Get coefficient estimates and standard errors
  coef_values <- model_summary$coefficients[, 1]  # Coefficients
  coef_se <- model_summary$coefficients[, 2]  # Standard errors
  
  # Determine the indices of selected variables
  term_labels <- attr(model$terms, "term.labels")
  ind <- match(term_labels, colnames(data)[-1])
  
  # Reorder coefficients and standard errors to match variable indices
  reordered_coef <- numeric(ncol(data))
  reordered_se <- numeric(ncol(data))
  
  reordered_coef[c(1, ind + 1)] <- coef_values
  reordered_se[c(1, ind + 1)] <- coef_se
  
  # Replace zero coefficients with NA
  reordered_coef[reordered_coef == 0] <- NA
  reordered_se[reordered_se == 0] <- NA
  
  # Create a named list to return results
  result <- list(
    coefficients = reordered_coef,
    standard_errors = reordered_se
  )
  
  return(result)
}

# Applying to above function to many different stability selection outputs
ApplySSModel <- function(Selection.Regs, Response, data) {
  results <- list()
  
  for (method in names(Selection.Regs)) {
    for (dataset in names(Selection.Regs[[method]])) {
      regressors <- Selection.Regs[[method]][[dataset]]
      results[[method]][[dataset]] <- FitStabilitySelectedModel(Response = Response, Regressors = regressors, data = data)
    }
  }
  
  return(results)
}
