### Functions and Packages
## Load packages
# If using cluster probably need to install `mvtnorm`
library(glmnet)
library(tidyverse)
library(patchwork)
library(parallel)
library(RColorBrewer)

# Function to create empty objects
Empty_Objects <- function(n, dims, Type = "matrix", ObjectNames = NULL, AssignGlobal = TRUE) {
  # Create list of all empty objects
  objects <- vector("list", n)
  
  # Create n objects based on which specified type
  for (i in seq_len(n)) {
    objects[[i]] <- switch(Type,
                           "vector" = rep(NA, dims[1]),
                           "matrix" = matrix(NA, nrow = dims[1], ncol = dims[2]),
                           "array" = matrix(NA, dim = c(dims[1], dims[2], dims[3])),
                           "data.frame" = as.data.frame(matrix(NA, nrow = dims[1], ncol = dims[2])),
                           stop("Unsupported type. Use 'vector', 'matrix', 'array', or 'data.frame'"))
  }
  # Name each element of list
  if (is.null(ObjectNames)) {
    names(objects) <- paste0(Type, "_", seq_len(n))
  } else {
    names(objects) <- ObjectNames
  }
  
  # Assign each list element as an individual object in the global environment
  if (AssignGlobal) {
    list2env(objects, envir = .GlobalEnv)
  }
  
  return(objects)
}

# Sequential Model Selection Function
ModelSelection.Fun <- function(mu.vec = rep(1, 3),          # Mean vector for covariates
                               cov.mat = diag(3),           # Covariance matrix
                               n = 200,                     # Sample Size
                               epsilon.error = 2,           # Error variance
                               beta.seq = 0:3,              # Coefficients (including intercept)
                               MS.type = "forward",         # Model selection type (forward, backwards, step-wise)
                               nsim = 1000,                 # Number of simulations
                               criteria = "AIC",            # Model selection criterion (AIC or BIC)
                               iterseed = TRUE,             # Set a new seed each iteration for easier traceback
                               verbose = FALSE){            # Print progress if desired
  ### Empty object creation
  Reg.Sel.List <- rep(NA, nsim)
  
  # Objects for model selected model
  res.SE <- rep(NA, nsim)
  Coef.mat <- CoefSE.mat <- Tval.mat <- Pval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CI.array <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  # Objects for "correct" model
  Cres.SE <- rep(NA, nsim)
  CCoef.mat <- CCoefSE.mat <- CTval.mat <- CPval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CCI.array <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  # Objects for full model
  Fres.SE <- rep(NA, nsim)
  FCoef.mat <- FCoefSE.mat <- FTval.mat <- FPval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  FCI.array <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  ### Start Simulations
  for(i in 1:nsim){
    ### Data Generation
    if (iterseed) {
      set.seed(i)
    }
    
    ## Random MVN Sample
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    
    # ABC labeling of regressors unless there is a lot
    if(length(mu.vec) > 22){
      colnames(x.mat) <- paste("X", 1:length(mu.vec), sep = "")
    }
    colnames(x.mat) <- LETTERS[1:length(mu.vec)]
    
    # Creating response
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, mean = 0, sd = sqrt(epsilon.error))
    
    # Creating dataframe
    da.data <- cbind.data.frame(Y = as.vector(y.vec), x.mat)
    
    # Perform model selection
    if(MS.type %in% c("forward", "backward", "both")) {
      ## Training Data
      # Null model for model selection (on training data)
      null.model <- lm(Y ~ 1, data = da.data)
      
      # Full model for model selection (on training data)
      full.model <- lm(Y ~ ., data = da.data)
      
      # Scope range
      Scope <- list(lower = as.formula("~ 1"),
                    upper = as.formula(paste("~", paste(colnames(x.mat), collapse = " + "))))
      
      # Depending on which model selection method
      if (MS.type == "forward") {
        Mod.Select <- step(object = null.model, scope = Scope[["upper"]], 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
        
      } else if (MS.type == "backward") {
        Mod.Select <- step(object = full.model, scope = Scope, 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
        
      } else {
        Mod.Select <- step(object = full.model, scope = Scope, 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
      }
      
      # Noting which regressors were chosen in the selected model
      Reg.Sel.List[i] <- paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = "")
      
      # Pulling summary from selected model
      lm.sel <- lm(formula = Mod.Select$call$formula, data = da.data)
      lm.select <- summary(lm.sel)
      
      # Checking which index the variables chosen correspond to
      ind <- match(attr(lm.select$terms, "term.labels"), colnames(x.mat))
      
      # Saving coefficients, coefficient standard errors, t-statistics, and model's residual SE
      Coef.mat[i, c(1, ind + 1)] <- lm.select$coefficients[, 1]
      CoefSE.mat[i, c(1, ind + 1)] <- lm.select$coefficients[, 2]
      Tval.mat[i, c(1, ind + 1)] <- lm.select$coefficients[, 3]
      Pval.mat[i, c(1, ind + 1)] <- lm.select$coefficients[, 4]
      CI.array[c(1, ind + 1), , i] <- confint(lm.sel)
      res.SE[i] <- lm.select$sigma
      
      # Pulling summary from "correct" model
      Correct.Regs <- colnames(x.mat)[which(beta.seq[-1] != 0)]
      lm.cor <- lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data)
      lm.correct <- summary(lm.cor)
      ind <- match(attr(lm.correct$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the correct model
      CCoef.mat[i, c(1, ind + 1)] <- lm.correct$coefficients[, 1]
      CCoefSE.mat[i, c(1, ind + 1)] <- lm.correct$coefficients[, 2]
      CTval.mat[i, c(1, ind + 1)] <- lm.correct$coefficients[, 3]
      CPval.mat[i, c(1, ind + 1)] <- lm.correct$coefficients[, 4]
      CCI.array[c(1, ind + 1), , i] <- confint(lm.cor)
      Cres.SE[i] <- lm.correct$sigma
      
      # Pulling summary from full model
      lm.fu <- lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = da.data)
      lm.full <- summary(lm.fu)
      ind <- match(attr(lm.full$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the full model
      FCoef.mat[i, c(1, ind + 1)] <- lm.full$coefficients[, 1]
      FCoefSE.mat[i, c(1, ind + 1)] <- lm.full$coefficients[, 2]
      FTval.mat[i, c(1, ind + 1)] <- lm.full$coefficients[, 3]
      FPval.mat[i, c(1, ind + 1)] <- lm.full$coefficients[, 4]
      FCI.array[c(1, ind + 1), , i] <- confint(lm.fu)
      Fres.SE[i] <- lm.full$sigma
    }
    
    # Print progress if desired
    if (verbose) {
      if(i %% (nsim / 10) == 0){
        cat(i, "Simulations completed", "\n")
      }
    }
  }
  
  ### Aggregating Results
  # Combine results for the model selected
  Results <- cbind.data.frame(Reg.Sel.List, Coef.mat, CoefSE.mat, Tval.mat, Pval.mat, res.SE)
  
  # Combine results for the "correct" model
  CResults <- cbind.data.frame(rep(paste(colnames(x.mat)[which(beta.seq[-1] != 0)], collapse = ""), nsim),
                               CCoef.mat, CCoefSE.mat, CTval.mat, CPval.mat, Cres.SE)
  
  # Combine results for the full model
  FResults <- cbind.data.frame(rep(paste(colnames(x.mat), collapse = ""), nsim),
                               FCoef.mat, FCoefSE.mat, FTval.mat, FPval.mat, Fres.SE)
  
  # Labeling 
  colnames(Results) <- c("Regressors", 
                         paste("Beta", 0:(length(beta.seq)-1), sep = ""), 
                         paste("SE.Beta", 0:(length(beta.seq)-1), sep = ""),
                         paste("T", 0:(length(beta.seq)-1), sep = ""),
                         paste("p.val", 0:(length(beta.seq) - 1), sep = ""),
                         "Res.SE")
  
  colnames(FResults) <- colnames(CResults) <- colnames(Results)
  
  # Labeling confidence interval arrays
  dimnames(CI.array) <- dimnames(CCI.array) <- dimnames(FCI.array) <- 
    list(dimnames(confint(lm.fu))[[1]],
         dimnames(confint(lm.fu))[[2]],
         paste("Sim", 1:nsim, sep = ""))
  
  
  # Output
  return(Output = list(Post.Model.Selection = Results,
                       Post.Model.Selection.CIs = CI.array,
                       Correct.Model = CResults, 
                       Correct.CIs = CCI.array,
                       Full.Model = FResults, 
                       Full.CIs = FCI.array))
}

# Parallel of above
ModelSelection.Parallel.Fun <- function(mu.vec = rep(1, 3),
                                        cov.mat = diag(3),
                                        n = 200,
                                        epsilon.error = 2,
                                        beta.seq = 0:3,
                                        MS.type = "forward",
                                        nsim = 1000,
                                        criteria = "AIC",
                                        iterseed = TRUE,
                                        iseed = 123,
                                        verbose = FALSE,
                                        ncores = parallel::detectCores() - 1) {
  
  # Prepare cluster
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  
  # Set reproducible RNG streams BEFORE running simulations
  if (iterseed) parallel::clusterSetRNGStream(cl, iseed = iseed)
  
  # Export variables and load required packages on each worker
  parallel::clusterExport(cl, varlist = c("mu.vec", "cov.mat", "n", "epsilon.error", "beta.seq",
                                          "MS.type", "criteria"), envir = environment())
  parallel::clusterEvalQ(cl, library(mvtnorm))
  
  run_simulation <- function(i) {
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    colnames(x.mat) <- if (length(mu.vec) > 22) paste0("X", seq_along(mu.vec)) else LETTERS[seq_along(mu.vec)]
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, sd = sqrt(epsilon.error))
    da.data <- data.frame(Y = as.vector(y.vec), x.mat)
    
    null.model <- lm(Y ~ 1, data = da.data)
    full.model <- lm(Y ~ ., data = da.data)
    Scope <- list(lower = ~1, upper = as.formula(paste("~", paste(colnames(x.mat), collapse = "+"))))
    
    k_val <- ifelse(criteria == "AIC", 2, log(n))
    Mod.Select <- switch(MS.type,
                         forward = step(null.model, scope = Scope$upper, direction = "forward", trace = 0, k = k_val),
                         backward = step(full.model, scope = Scope, direction = "backward", trace = 0, k = k_val),
                         both = step(full.model, scope = Scope, direction = "both", trace = 0, k = k_val))
    
    lm.sel <- lm(formula = Mod.Select$call$formula, data = da.data)
    lm.select <- summary(lm.sel)
    ind.sel <- match(attr(lm.select$terms, "term.labels"), colnames(x.mat))
    
    Correct.Regs <- colnames(x.mat)[which(beta.seq[-1] != 0)]
    lm.cor <- lm(as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data)
    lm.correct <- summary(lm.cor)
    ind.cor <- match(attr(lm.correct$terms, "term.labels"), colnames(x.mat))
    
    lm.fu <- lm(as.formula(paste("Y ~", paste(colnames(x.mat), collapse = "+"))), data = da.data)
    lm.full <- summary(lm.fu)
    ind.full <- match(attr(lm.full$terms, "term.labels"), colnames(x.mat))
    
    CI.array <- array(NA, dim = c(length(beta.seq), 2))
    CI.array[c(1, ind.sel + 1), ] <- confint(lm.sel)
    
    CCI.array <- array(NA, dim = c(length(beta.seq), 2))
    CCI.array[c(1, ind.cor + 1), ] <- confint(lm.cor)
    
    FCI.array <- array(NA, dim = c(length(beta.seq), 2))
    FCI.array[c(1, ind.full + 1), ] <- confint(lm.fu)
    
    list(
      RegSel = paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = ""),
      Coef = replace(rep(NA, length(beta.seq)), c(1, ind.sel + 1), lm.select$coefficients[, 1]),
      SE = replace(rep(NA, length(beta.seq)), c(1, ind.sel + 1), lm.select$coefficients[, 2]),
      Tval = replace(rep(NA, length(beta.seq)), c(1, ind.sel + 1), lm.select$coefficients[, 3]),
      Pval = replace(rep(NA, length(beta.seq)), c(1, ind.sel + 1), lm.select$coefficients[, 4]),
      ResSE = lm.select$sigma,
      CI = CI.array,
      
      CCoef = replace(rep(NA, length(beta.seq)), c(1, ind.cor + 1), lm.correct$coefficients[, 1]),
      CSE = replace(rep(NA, length(beta.seq)), c(1, ind.cor + 1), lm.correct$coefficients[, 2]),
      CTval = replace(rep(NA, length(beta.seq)), c(1, ind.cor + 1), lm.correct$coefficients[, 3]),
      CPval = replace(rep(NA, length(beta.seq)), c(1, ind.cor + 1), lm.correct$coefficients[, 4]),
      CresSE = lm.correct$sigma,
      CCI = CCI.array,
      
      FCoef = replace(rep(NA, length(beta.seq)), c(1, ind.full + 1), lm.full$coefficients[, 1]),
      FSE = replace(rep(NA, length(beta.seq)), c(1, ind.full + 1), lm.full$coefficients[, 2]),
      FTval = replace(rep(NA, length(beta.seq)), c(1, ind.full + 1), lm.full$coefficients[, 3]),
      FPval = replace(rep(NA, length(beta.seq)), c(1, ind.full + 1), lm.full$coefficients[, 4]),
      FresSE = lm.full$sigma,
      FCI = FCI.array
    )
  }
  
  sim_results <- parallel::parLapply(cl, 1:nsim, run_simulation)
  
  # Aggregation
  Reg.Sel.List <- sapply(sim_results, `[[`, "RegSel")
  Coef.mat     <- do.call(rbind, lapply(sim_results, `[[`, "Coef"))
  CoefSE.mat   <- do.call(rbind, lapply(sim_results, `[[`, "SE"))
  Tval.mat     <- do.call(rbind, lapply(sim_results, `[[`, "Tval"))
  Pval.mat     <- do.call(rbind, lapply(sim_results, `[[`, "Pval"))
  res.SE       <- sapply(sim_results, `[[`, "ResSE")
  
  CCoef.mat    <- do.call(rbind, lapply(sim_results, `[[`, "CCoef"))
  CCoefSE.mat  <- do.call(rbind, lapply(sim_results, `[[`, "CSE"))
  CTval.mat    <- do.call(rbind, lapply(sim_results, `[[`, "CTval"))
  CPval.mat    <- do.call(rbind, lapply(sim_results, `[[`, "CPval"))
  Cres.SE      <- sapply(sim_results, `[[`, "CresSE")
  
  FCoef.mat    <- do.call(rbind, lapply(sim_results, `[[`, "FCoef"))
  FCoefSE.mat  <- do.call(rbind, lapply(sim_results, `[[`, "FSE"))
  FTval.mat    <- do.call(rbind, lapply(sim_results, `[[`, "FTval"))
  FPval.mat    <- do.call(rbind, lapply(sim_results, `[[`, "FPval"))
  Fres.SE      <- sapply(sim_results, `[[`, "FresSE")
  
  CI.array <- simplify2array(lapply(sim_results, `[[`, "CI"))
  CCI.array <- simplify2array(lapply(sim_results, `[[`, "CCI"))
  FCI.array <- simplify2array(lapply(sim_results, `[[`, "FCI"))
  
  Results <- cbind.data.frame(Reg.Sel.List, Coef.mat, CoefSE.mat, Tval.mat, Pval.mat, res.SE)
  CResults <- cbind.data.frame(rep(paste(colnames(mu.vec), collapse = ""), nsim),
                               CCoef.mat, CCoefSE.mat, CTval.mat, CPval.mat, Cres.SE)
  FResults <- cbind.data.frame(rep(paste(colnames(mu.vec), collapse = ""), nsim),
                               FCoef.mat, FCoefSE.mat, FTval.mat, FPval.mat, Fres.SE)
  
  colnames(Results) <- colnames(FResults) <- colnames(CResults) <- c(
    "Regressors",
    paste0("Beta", 0:(length(beta.seq)-1)),
    paste0("SE.Beta", 0:(length(beta.seq)-1)),
    paste0("T", 0:(length(beta.seq)-1)),
    paste0("p.val", 0:(length(beta.seq)-1)),
    "Res.SE"
  )
  
  dimnames(CI.array) <- dimnames(CCI.array) <- dimnames(FCI.array) <-
    list(paste0("Beta", 0:(length(beta.seq)-1)), c("2.5 %", "97.5 %"), paste0("Sim", 1:nsim))
  
  return(list(Post.Model.Selection = Results,
              Post.Model.Selection.CIs = CI.array,
              Correct.Model = CResults,
              Correct.CIs = CCI.array,
              Full.Model = FResults,
              Full.CIs = FCI.array))
}

# Same function as above but not parallel and with data splitting
ModelSelection.DataSplit.Fun <- function(mu.vec = rep(0, 3),          # Mean vector for covariates
                                         cov.mat = diag(3),           # Covariance matrix
                                         n = 200,                     # Sample Size
                                         epsilon.error = 2,           # Error variance
                                         beta.seq = 0:3,              # Coefficients (including intercept)
                                         MS.type = "forward",         # Model selection type (forward, backwards, step-wise)
                                         nsim = 1000,                 # Number of simulations
                                         criteria = "AIC",            # Model selection criterion (AIC or BIC)
                                         SP = 0.70,                   # Data splitting ratio (training data proportion)
                                         iterseed = TRUE,             # Set a new seed each iteration for easier traceback
                                         verbose = FALSE){            # Print progress if desired
  ### Empty object creation
  ## Initialize objects for train and test data
  Reg.Sel.List.train <- rep(NA, nsim)
  Reg.Sel.List.test <- rep(NA, nsim)
  
  # Objects for model selected model for train and test
  res.SE.train <- rep(NA, nsim)
  res.SE.test <- rep(NA, nsim)
  
  Coef.mat.train <- CoefSE.mat.train <- Tval.mat.train <- Pval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CI.array.train <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  Coef.mat.test <- CoefSE.mat.test <- Tval.mat.test <- Pval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CI.array.test <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  # Objects for "correct" model for train and test
  Cres.SE.train <- rep(NA, nsim)
  Cres.SE.test <- rep(NA, nsim)
  
  CCoef.mat.train <- CCoefSE.mat.train <- CTval.mat.train <- CPval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CCI.array.train <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  CCoef.mat.test <- CCoefSE.mat.test <- CTval.mat.test <- CPval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  CCI.array.test <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  # Objects for full model for train and test
  Fres.SE.train <- rep(NA, nsim)
  Fres.SE.test <- rep(NA, nsim)
  
  FCoef.mat.train <- FCoefSE.mat.train <- FTval.mat.train <- FPval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  FCI.array.train <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  FCoef.mat.test <- FCoefSE.mat.test <- FTval.mat.test <- FPval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  FCI.array.test <- array(data = NA, dim = c(length(beta.seq), 2, nsim))
  
  ### Start Simulations
  for(i in 1:nsim){
    ### Data Generation
    if (iterseed) {
      set.seed(i)
    }
    
    ## Random MVN Sample
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    
    # ABC labeling of regressors unless there is a lot
    if(length(mu.vec) > 22){
      colnames(x.mat) <- paste("X", 1:length(mu.vec), sep = "")
    }
    colnames(x.mat) <- LETTERS[1:length(mu.vec)]
    
    # Creating response
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, mean = 0, sd = sqrt(epsilon.error))
    
    # Creating dataframe
    da.data <- cbind.data.frame(Y = as.vector(y.vec), x.mat)
    
    # Data-Splitting Strategy - Performing specified model selection technique
    # Split data into training and test sets
    index <- sample(1:n, size = round(n * SP))
    train.data <- da.data[index, ]
    test.data <- da.data[-index, ]
    
    # Perform model selection only on training data
    if(MS.type %in% c("forward", "backward", "both")) {
      ## Training Data
      # Null model for model selection (on training data)
      null.model <- lm(Y ~ 1, data = train.data)
      
      # Full model for model selection (on training data)
      full.model <- lm(Y ~ ., data = train.data)
      
      # Scope range
      Scope <- list(lower = as.formula("~ 1"),
                    upper = as.formula(paste("~", paste(colnames(x.mat), collapse = " + "))))
      
      # Depending on which model selection method
      if (MS.type == "forward") {
        Mod.Select <- step(object = null.model, scope = Scope[["upper"]], 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
        
      } else if (MS.type == "backward") {
        Mod.Select <- step(object = full.model, scope = Scope, 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
        
      } else {
        Mod.Select <- step(object = full.model, scope = Scope, 
                           direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
      }
      
      # Noting which regressors were chosen in the selected model (for train data)
      Reg.Sel.List.train[i] <- paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = "")
      
      # Pulling summary from selected model (on train data)
      lm.sel <- lm(formula = Mod.Select$call$formula, data = train.data)
      lm.select <- summary(lm.sel)
      
      # Checking which index the variables chosen correspond to
      ind <- match(attr(lm.select$terms, "term.labels"), colnames(x.mat))
      
      # Saving coefficients, coefficient standard errors, t-statistics, and model's residual SE (for train data)
      Coef.mat.train[i, c(1, ind + 1)] <- lm.select$coefficients[, 1]
      CoefSE.mat.train[i, c(1, ind + 1)] <- lm.select$coefficients[, 2]
      Tval.mat.train[i, c(1, ind + 1)] <- lm.select$coefficients[, 3]
      Pval.mat.train[i, c(1, ind + 1)] <- lm.select$coefficients[, 4]
      CI.array.train[c(1, ind + 1), , i] <- confint(lm.sel)
      res.SE.train[i] <- lm.select$sigma
      
      # Pulling summary from "correct" model (on train data)
      Correct.Regs <- colnames(x.mat)[which(beta.seq[-1] != 0)]
      lm.cor <- lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = train.data)
      lm.correct <- summary(lm.cor)
      ind <- match(attr(lm.correct$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the correct model (for train data)
      CCoef.mat.train[i, c(1, ind + 1)] <- lm.correct$coefficients[, 1]
      CCoefSE.mat.train[i, c(1, ind + 1)] <- lm.correct$coefficients[, 2]
      CTval.mat.train[i, c(1, ind + 1)] <- lm.correct$coefficients[, 3]
      CPval.mat.train[i, c(1, ind + 1)] <- lm.correct$coefficients[, 4]
      CCI.array.train[c(1, ind + 1), , i] <- confint(lm.cor)
      Cres.SE.train[i] <- lm.correct$sigma
      
      # Pulling summary from full model (on train data)
      lm.fu <- lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = train.data)
      lm.full <- summary(lm.fu)
      ind <- match(attr(lm.full$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the full model (for train data)
      FCoef.mat.train[i, c(1, ind + 1)] <- lm.full$coefficients[, 1]
      FCoefSE.mat.train[i, c(1, ind + 1)] <- lm.full$coefficients[, 2]
      FTval.mat.train[i, c(1, ind + 1)] <- lm.full$coefficients[, 3]
      FPval.mat.train[i, c(1, ind + 1)] <- lm.full$coefficients[, 4]
      FCI.array.train[c(1, ind + 1), , i] <- confint(lm.fu)
      Fres.SE.train[i] <- lm.full$sigma
      
      ### Test Data
      ## Apply the selected model to the test data
      # Check if Mod.Select$terms has any terms
      selected_terms <- sort(attr(Mod.Select$terms, "term.labels"))
      if (length(selected_terms) > 0) {
        # Create the formula string
        formula_string <- paste("Y ~", paste(selected_terms, collapse = " + "))
        
        # Create the model
        test.model <- lm(as.formula(formula_string), data = test.data)
        
      } else {
        next
      }
      
      # Pulling summary from selected model (on test data)
      lm.select.test <- summary(test.model)
      
      # Checking which index the variables chosen correspond to
      ind.test <- match(attr(test.model$terms, "term.labels"), colnames(x.mat))
      
      # Saving coefficients, coefficient standard errors, t-statistics, and model's residual SE (for test data)
      Coef.mat.test[i, c(1, ind.test + 1)] <- lm.select.test$coefficients[, 1]
      CoefSE.mat.test[i, c(1, ind.test + 1)] <- lm.select.test$coefficients[, 2]
      Tval.mat.test[i, c(1, ind.test + 1)] <- lm.select.test$coefficients[, 3]
      Pval.mat.test[i, c(1, ind.test + 1)] <- lm.select.test$coefficients[, 4]
      CI.array.test[c(1, ind + 1), , i] <- confint(test.model)
      res.SE.test[i] <- lm.select.test$sigma
      
      # Pulling summary from "correct" model (on test data)
      lm.cor.test <- lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = test.data)
      lm.correct.test <- summary(lm.cor.test)
      ind.test <- match(attr(lm.correct.test$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the correct model (for test data)
      CCoef.mat.test[i, c(1, ind.test + 1)] <- lm.correct.test$coefficients[, 1]
      CCoefSE.mat.test[i, c(1, ind.test + 1)] <- lm.correct.test$coefficients[, 2]
      CTval.mat.test[i, c(1, ind.test + 1)] <- lm.correct.test$coefficients[, 3]
      CPval.mat.test[i, c(1, ind.test + 1)] <- lm.correct.test$coefficients[, 4]
      CCI.array.test[c(1, ind + 1), , i] <- confint(lm.cor.test)
      Cres.SE.test[i] <- lm.correct.test$sigma
      
      # Pulling summary from full model (on test data)
      lm.fu.test <- lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = test.data)
      lm.full.test <- summary(lm.fu.test)
      ind.test <- match(attr(lm.full.test$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the full model (for test data)
      FCoef.mat.test[i, c(1, ind.test + 1)] <- lm.full.test$coefficients[, 1]
      FCoefSE.mat.test[i, c(1, ind.test + 1)] <- lm.full.test$coefficients[, 2]
      FTval.mat.test[i, c(1, ind.test + 1)] <- lm.full.test$coefficients[, 3]
      FPval.mat.test[i, c(1, ind.test + 1)] <- lm.full.test$coefficients[, 4]
      FCI.array.test[c(1, ind + 1), , i] <- confint(lm.fu.test)
      Fres.SE.test[i] <- lm.full.test$sigma
    }
    
    # Print progress if desired
    if (verbose) {
      if(i %% (nsim / 10) == 0){
        cat(i, "Simulations completed", "\n")
      }
    }
  }
  
  ### Aggregating Results
  # Combine results for the model selected from the training data
  Train.Results <- cbind.data.frame(Reg.Sel.List.train, Coef.mat.train, CoefSE.mat.train, Tval.mat.train, Pval.mat.train, res.SE.train)
  Test.Results <- cbind.data.frame(Reg.Sel.List.test, Coef.mat.test, CoefSE.mat.test, Tval.mat.test, Pval.mat.test, res.SE.test)
  
  # Combine results for the "correct" model
  Train.CResults <- cbind.data.frame(CCoef.mat.train, CCoefSE.mat.train, CTval.mat.train, CPval.mat.train, Cres.SE.train)
  Test.CResults <- cbind.data.frame(CCoef.mat.test, CCoefSE.mat.test, CTval.mat.test, CPval.mat.test, Cres.SE.test)
  
  # Combine results for the full model
  Train.FResults <- cbind.data.frame(FCoef.mat.train, FCoefSE.mat.train, FTval.mat.train, FPval.mat.train, Fres.SE.train)
  Test.FResults <- cbind.data.frame(FCoef.mat.test, FCoefSE.mat.test, FTval.mat.test, FPval.mat.test, Fres.SE.test)
  
  # Labeling for Train and Test Results
  colnames(Train.Results) <- c("Regressors", paste("Beta", 0:(length(beta.seq) - 1), sep = ""), 
                               paste("SE.Beta", 0:(length(beta.seq) - 1), sep = ""),
                               paste("T.val", 0:(length(beta.seq) - 1), sep = ""),
                               paste("p.val", 0:(length(beta.seq) - 1), sep = ""),
                               "Res.SE")
  colnames(Test.Results) <- colnames(Train.Results)
  
  # Labeling confidence interval arrays
  dimnames(CI.array.train) <- dimnames(CI.array.test) <- dimnames(CCI.array.train) <- dimnames(CCI.array.test) <-
    dimnames(FCI.array.train) <- dimnames(FCI.array.test) <- list(dimnames(confint(lm.fu.test))[[1]],
                                                                  dimnames(confint(lm.fu.test))[[2]],
                                                                  paste("Sim", 1:nsim, sep = ""))
  
  # Labeling for Correct Model Results (Train and Test)
  colnames(Train.CResults) <- colnames(Train.Results)[-1]
  colnames(Test.CResults) <- colnames(Train.Results)[-1]
  
  # Labeling for Full Model Results (Train and Test)
  colnames(Train.FResults) <- colnames(Train.Results)[-1]
  colnames(Test.FResults) <- colnames(Train.Results)[-1]
  
  # Output
  return(Output = list(
    Train.Post.Model.Selection = Train.Results,
    Train.Post.Model.Selection.CIs = CI.array.train,
    Test.Post.Model.Selection = Test.Results,
    Test.Post.Model.Selection.CIs = CI.array.test,
    Train.Correct.Model = Train.CResults,
    Train.Correct.CIs = CCI.array.train,
    Test.Correct.Model = Test.CResults,
    Test.Correct.CIs = CCI.array.test,
    Train.Full.Model = Train.FResults,
    Train.Full.CIs = FCI.array.train,
    Test.Full.Model = Test.FResults,
    Test.Full.CIs = FCI.array.test))
  
}

# LASSO and Ridge Model Selection Function
Penalized.MS.Fun <- function(mu.vec = rep(0, 3),          # Mean vector for covariates
                             cov.mat = diag(3),           # Covariance matrix
                             n = 200,                     # Sample Size
                             epsilon.error = 2,           # Error variance
                             beta.seq = 0:3,              # Coefficients (including intercept)
                             nsim = 1000,                 # Number of simulations
                             iterseed = TRUE,             # Set a new seed each iteration for easier traceback
                             verbose = FALSE){            # Print progress if desired
  # Objects for LASSO and Ridge Models
  LCoef.mat <- RCoef.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
  
  # Objects for LASSO and Ridge Models for train and test
  # LCoef.mat.train <- RCoef.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
  # LCoef.mat.test <- RCoef.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
  
  ### Start Simulations
  for(i in 1:nsim){
    ### Data Generation
    if (iterseed) {
      set.seed(i)
    }
    
    ## Random MVN Sample
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    
    # ABC labeling of regressors unless there is a lot
    if(length(mu.vec) > 22){
      colnames(x.mat) <- paste("X", 1:length(mu.vec), sep = "")
    }
    colnames(x.mat) <- LETTERS[1:length(mu.vec)]
    
    # Creating response
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, mean = 0, sd = sqrt(epsilon.error))
    
    # Creating dataframe
    da.data <- cbind.data.frame(Y = as.vector(y.vec), x.mat)
    
    ### Penalized Regression
    # LASSO 
    # Cross validation to calculate lambda
    LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    
    # Fitting Model
    LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    LCoef.mat[i, 1:(length(beta.seq) - 1)] <- as.matrix(LASSO.mod$beta)
    
    # Ridge
    # Cross validation to calculate lambda
    Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    
    # Fitting Model
    Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                        standardize = TRUE, thresh = 1e-10)
    # Pulling coefficients
    RCoef.mat[i, 1:(length(beta.seq) - 1)] <- as.matrix(Ridge.mod$beta)
    
    # Print progress if desired
    if (verbose) {
      if(i %% (nsim / 10) == 0){
        cat(i, "Simulations completed", "\n")
      }
    }
  }
  
  # Adding a regressor column
  LRegs <- apply(LCoef.mat, 1, function(x) paste(colnames(x.mat)[which(x != 0)], collapse = ""))
  L.Coef.df <- data.frame(Regressors = LRegs, LCoef.mat)
  
  R.Coef.df <- data.frame(Regressors = rep(paste(colnames(x.mat), collapse = ""), nsim),
                          RCoef.mat)
  
  # Labeling 
  colnames(L.Coef.df) <- colnames(R.Coef.df) <- c("Regressors", paste("Beta", 1:(length(beta.seq) - 1), sep = ""))
  
  # Output
  return(Output = list(LASSO = L.Coef.df,
                       Ridge = R.Coef.df))
}
# Parallel of above
Penalized.MS.Parallel.Fun <- function(mu.vec = rep(0, 3),
                                      cov.mat = diag(3),
                                      n = 200,
                                      epsilon.error = 2,
                                      beta.seq = 0:3,
                                      nsim = 1000,
                                      iterseed = TRUE,
                                      iseed = 123,
                                      verbose = FALSE,
                                      ncores = parallel::detectCores() - 1) {
  
  # Prepare cluster
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  
  # Set RNG stream
  if (iterseed) parallel::clusterSetRNGStream(cl, iseed = iseed)
  
  # Export necessary objects and libraries
  parallel::clusterExport(cl, varlist = c("mu.vec", "cov.mat", "n", "epsilon.error", "beta.seq"), envir = environment())
  parallel::clusterEvalQ(cl, { library(glmnet); library(mvtnorm) })
  
  # Simulation function
  run_simulation <- function(i) {
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    colnames(x.mat) <- if (length(mu.vec) > 22) paste0("X", seq_along(mu.vec)) else LETTERS[seq_along(mu.vec)]
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, mean = 0, sd = sqrt(epsilon.error))
    
    # LASSO
    LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min,
                        standardize = TRUE, thresh = 1e-10)
    LCoef <- as.numeric(LASSO.mod$beta)
    
    # Ridge
    Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
    Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min,
                        standardize = TRUE, thresh = 1e-10)
    RCoef <- as.numeric(Ridge.mod$beta)
    
    # Get variable names
    colnames_used <- colnames(x.mat)
    
    list(Lasso.Coef = LCoef, Ridge.Coef = RCoef, colnames = colnames_used)
  }
  
  # Run simulations
  sim_results <- parallel::parLapply(cl, 1:nsim, run_simulation)
  
  # Extract matrices
  LCoef.mat <- do.call(rbind, lapply(sim_results, function(x) x$Lasso.Coef))
  RCoef.mat <- do.call(rbind, lapply(sim_results, function(x) x$Ridge.Coef))
  varnames <- sim_results[[1]]$colnames
  
  # Construct regressor strings
  LRegs <- apply(LCoef.mat, 1, function(coefs) paste(varnames[which(coefs != 0)], collapse = ""))
  L.Coef.df <- data.frame(Regressors = LRegs, LCoef.mat)
  R.Coef.df <- data.frame(Regressors = rep(paste(varnames, collapse = ""), nsim), RCoef.mat)
  
  # Labeling
  colnames(L.Coef.df) <- colnames(R.Coef.df) <- c("Regressors", paste0("Beta", 1:(length(beta.seq) - 1)))
  
  return(list(LASSO = L.Coef.df, Ridge = R.Coef.df))
}
# Fits all possible models
All.Models.Parallel.Fun <- function(mu.vec = rep(1, 3),
                                    cov.mat = diag(3),
                                    n = 200,
                                    epsilon.error = 2,
                                    beta.seq = 0:3,
                                    nsim = 1000,
                                    iterseed = TRUE,
                                    iseed = 123,
                                    verbose = FALSE,
                                    ncores = parallel::detectCores() - 1) {
  
  # Prepare cluster
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  if (iterseed) parallel::clusterSetRNGStream(cl, iseed = iseed)
  
  # Export variables and load required packages
  parallel::clusterExport(cl, varlist = c("mu.vec", "cov.mat", "n", "epsilon.error", "beta.seq"), envir = environment())
  parallel::clusterEvalQ(cl, library(mvtnorm))
  
  run_simulation <- function(i) {
    x.mat <- mvtnorm::rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
    colnames(x.mat) <- if (length(mu.vec) > 22) paste0("X", seq_along(mu.vec)) else LETTERS[seq_along(mu.vec)]
    y.vec <- beta.seq[1] + x.mat %*% beta.seq[-1] + rnorm(n, sd = sqrt(epsilon.error))
    da.data <- data.frame(Y = as.vector(y.vec), x.mat)
    
    # Generate all combinations of the variables
    variables <- colnames(x.mat)
    combinations <- unlist(lapply(0:length(variables), function(k) combn(variables, k, simplify = FALSE)), recursive = FALSE)
    combinations[[1]] <- character(0) # Correct intercept-only model
    
    # Fit all models and store results
    model_results <- lapply(combinations, function(vars) {
      regressor_label <- paste(sort(vars), collapse = "")
      formula <- as.formula(paste("Y ~", if (length(vars) == 0) "1" else paste(vars, collapse = " + ")))
      lm.fit <- lm(formula, data = da.data)
      coef_values <- summary(lm.fit)$coefficients[, 1]
      
      # Create a full-length coefficient vector with NA
      if (length(vars) == 0) {
        all_coefs <- setNames(rep(NA, length(variables) + 1), c("(Intercept)", variables))
        all_coefs["(Intercept)"] <- coef_values
      } else {
        all_coefs <- setNames(rep(NA, length(variables) + 1), c("(Intercept)", variables))
        all_coefs[names(coef_values)] <- coef_values
      }
      
      data.frame(Regressor = regressor_label, t(all_coefs))
    })
    do.call(rbind, model_results)
  }
  
  sim_results <- parallel::parLapply(cl, 1:nsim, run_simulation)
  final_results <- do.call(rbind, sim_results)
  
  return(final_results)
}


# Post-Model Selection Coefficient Plotting
Plot.Beta.Fun <- function(SimResults, MS.Type = c("Forward", "Backward", "Both", "LASSO", "Ridge"), 
                          SampleSize = c("LowSS", "HighSS"), Error = c("LowError", "HighError"), 
                          Correlation = c("Independent", "ModerateCorrelation"), IC = c("AIC", "BIC"),
                          BetaSetting = c("PatternA", "PatternB"), plot.stat = "Beta1") {
  ### Default Settings
  MS.Type <- match.arg(MS.Type)
  SampleSize <- match.arg(SampleSize)
  Error <- match.arg(Error)
  Correlation <- match.arg(Correlation)
  IC <- match.arg(IC)
  BetaSetting <- match.arg(BetaSetting)
  
  # Penalized or Sequential Model Selection
  if (MS.Type == "LASSO" || MS.Type == "Ridge") {
    PlotData <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["LASSO.Ridge"]][[MS.Type]]
    
  } else {
    PlotData <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Post.Model.Selection"]]
  }
  
  # Load in all the fitted model coefficients
  AllMods <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["AllModels"]]
  colnames(AllMods) <- c("Regressors", paste("Beta", 0:(dim(AllMods)[2] - 2), sep = ""))

  # Finding ranges for x-axis
  Min <- sapply(select(PlotData, matches(paste0("^", plot.stat))), min, na.rm = TRUE)
  Max <- sapply(select(PlotData, matches(paste0("^", plot.stat))), max, na.rm = TRUE)
  
  # Extract unique covariates from the Regressors column
  Mods <- unique(PlotData$Regressors)
  
  # Subset the data for each covariate
  subsets_all <- subsets <- list()
  for (mod in Mods) {
    subsets[[mod]] <- subset(PlotData, Regressors == mod)
    subsets_all[[mod]] <- subset(AllMods, Regressors == mod)
  }
  
  if(MS.Type == "LASSO" || MS.Type == "Ridge") {
    # Turn 0's into NA's for plotting
    subsets <- lapply(subsets, function(x) {
      x[,-1][x[,-1] == 0] <- NA
      return(x)
    })
  }
  
  # Generate a color palette for the `fill`
  num_colors <- length(Mods)
  fill_colors <- c(brewer.pal(min(max(num_colors, 3), 9), "Set1"), brewer.pal(min(max(num_colors, 3), 8), "Set2"))

  # Create individual plots for each subset
  plot_index <- 0
  individual_plots <- list()
  for (i in seq_along(Mods)) {
    subset_data <- subsets[[Mods[i]]]
    subset_all_data <- subsets_all[[Mods[i]]]
    
    # Check if there is enough models to even plot
    if(sum(!is.na(subset_data[[plot.stat]])) < 5){
      next
    }
    plot_index <- plot_index + 1
    
    # Create the plot
    plot <- ggplot() +
      geom_density(data = subset_data, aes_string(x = plot.stat, fill = "Regressors"), alpha = 0.33, color = "black") +
      geom_density(data = subset_all_data, aes_string(x = plot.stat), color = "black", fill = NA) +
      scale_fill_manual(values = fill_colors[plot_index], name = "Regressor") +
      labs(title = paste("Model", Mods[i]),
           x = "",
           y = "Density") +
      scale_x_continuous(limits = c(-max(abs(Min), abs(Max))-0.1, max(abs(Min), abs(Max))+0.1)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )

    # Add vertical line for Beta coefficient if Beta.seq is provided
    if (!is.null(Beta) && startsWith(plot.stat, "Beta")) {
      beta_number <- as.numeric(sub("Beta", "", plot.stat))
      if (!is.na(beta_number) && beta_number + 1 <= length(Beta)) {
        plot <- plot +
          geom_vline(xintercept = Beta[beta_number + 1], linetype = "dashed", color = "black")
      } else {
        warning("Invalid Beta coefficient index or Beta.seq is too short.")
      }
    }
    
    # Add the number of observations in the top-right corner
    n_obs <- sum(!is.na(subset_data[[plot.stat]]))
    plot <- plot +
      annotate("text", x = Inf, y = Inf, label = paste("n =", n_obs), 
               hjust = 1.1, vjust = 1.1, size = 4, fontface = "bold")
    
    individual_plots[[length(individual_plots) + 1]] <- plot
  }
  
  # Trying to determine number of columns dynamically based on the number of plots
  n_plots <- length(individual_plots)
  
  if (n_plots == 0) {
    return("Not enough data for any plots")
    # return(subsets)
  }
  
  ncol <- switch(n_plots, 
                 1, 2, 3, 2, 3, 
                 3, 2, 2, 3, 2, 
                 3, 3, 4, 4, 4, 
                 4)
  
  # Add caption with the simulation settings
  if (MS.Type == "LASSO" || MS.Type == "Ridge") {
    caption_text <- paste("Simulation Settings:",
                          "Sample Size:", SampleSize, "|",
                          "Error Level:", Error, "|",
                          "Correlation Setting:", Correlation, "|",
                          "Model Selection Type:", MS.Type) 
  } else {
    caption_text <- paste("Simulation Settings:",
                          "Sample Size:", SampleSize, "|",
                          "Error Level:", Error, "|",
                          "Correlation Setting:", Correlation, "|",
                          "Information Criterion:", IC, "|",
                          "Model Selection Type:", MS.Type)
  }
  
  
  # Combine all individual plots using patchwork and add the caption
  combined_plot <- wrap_plots(individual_plots) + 
    plot_layout(ncol = ncol, guides = "collect") &
    plot_annotation(title = bquote(beta[.(beta_number)] ~ "Coefficient Estimates"),
                    caption = caption_text, 
                    theme = theme(plot.caption = element_text(size = 10, face = "italic", hjust = 0.5),
                                  plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
  
  # Return the combined plot
  print(combined_plot)
  return(combined_plot)
}

# Don't really want to call it 'n' in the plots because it's not really the sample size
# Not currently compatible with data-splitting output.
# Perhaps overlay the full or correct model in each plot?
# I'd like each color to correspond to a specific plot / model...
# Some plot's tails are still getting cut off....????

# For odd numbers like 3, 5, 7, etc can use layout to make it look not terrible
Plot.Beta.Mix.Fun <- function(SimResults, MS.Type = c("Forward", "Backward", "Both"), 
                              SampleSize = c("LowSS", "HighSS"), Error = c("LowError", "HighError"), 
                              Correlation = c("Independent", "ModerateCorrelation"), IC = c("AIC", "BIC"),
                              BetaSetting = c("PatternA", "PatternB"), LR = TRUE, Beta.seq = NULL) {
  ### Default Settings
  MS.Type <- match.arg(MS.Type)
  SampleSize <- match.arg(SampleSize)
  Error <- match.arg(Error)
  Correlation <- match.arg(Correlation)
  IC <- match.arg(IC)
  BetaSetting <- match.arg(BetaSetting)
  
  # Data Extraction
  if (LR) {
    PlotDataL <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["LASSO.Ridge"]][["LASSO"]]
    PlotDataR <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["LASSO.Ridge"]][["Ridge"]]
    PlotDataPost <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Post.Model.Selection"]]
    PlotDataFull <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Full.Model"]]
    PlotDataCorrect <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Correct.Model"]]
    
  } else {
    PlotDataPost <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Post.Model.Selection"]]
    PlotDataFull <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Full.Model"]]
    PlotDataCorrect <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Correct.Model"]]
  }
  
  # Beta to plot (skipping Beta0)
  if (LR) {
    PlotBetas <- paste("Beta", 1:(length(Beta.seq) - 1), sep = "")
  } else {
    PlotBetas <- paste("Beta", 0:(length(Beta.seq) - 1), sep = "")
  }
  
  # Create individual plots for each beta
  individual_plots <- list()
  for(beta in PlotBetas){
    PlottingData <- data.frame(BetaPost = PlotDataPost[[beta]],
                               BetaFull = PlotDataFull[[beta]],
                               BetaCorr = PlotDataCorrect[[beta]])
    
    if (LR) {
      PlottingData <- data.frame(BetaPost = PlotDataPost[[beta]],
                                 BetaFull = PlotDataFull[[beta]],
                                 BetaCorr = PlotDataCorrect[[beta]],
                                 BetaLASSO = PlotDataL[[beta]],
                                 BetaRidge = PlotDataR[[beta]])
    }
    
    plot_title <- bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", beta)))] ~ " Sampling Distribution")
    plot_x_label <- bquote(hat(beta)[.(as.numeric(sub("Beta", "", beta)))])
    
    # Create the mixture density plot
    plot <- ggplot(PlottingData) +
      geom_density(aes(x = BetaPost, fill = "Selected Model"), alpha = 0.333, color = "#5ab9e8") +
      geom_density(aes(x = BetaFull, fill = "Full Model"), alpha = 0.333, color = "#7beb78") +
      geom_density(aes(x = BetaCorr, fill = "Correct Model"), alpha = 0.333, color = "#d55ae8") +
      {
        if (LR) {
          geom_density(aes(x = BetaLASSO, fill = "LASSO Model"), alpha = 0.333, color = "#fc031c")
        }
      } +
      {
        if (LR) {
          geom_density(aes(x = BetaRidge, fill = "Ridge Model"), alpha = 0.333, color = "#ff983d")
        }
      }
      
    
    # Count the number of non-NA columns (corresponding to present legend elements)
    num_legends <- sum(colSums(!is.na(PlottingData)) > 0)

    # Conditionally add the legend based on the number of unique elements
    if (LR) {
      if(num_legends >= 5){
        plot <- plot + 
          scale_fill_manual(name = "", values = c("Selected Model" = "#5ab9e8", 
                                                  "Full Model" = "#7beb78", 
                                                  "Correct Model" = "#d55ae8",
                                                  "LASSO Model" = "#fc031c",
                                                  "Ridge Model" = "#ff983d"))
      } else {
        plot <- plot + 
          scale_fill_manual(name = "", values = c("Selected Model" = "#5ab9e8", 
                                                  "Full Model" = "#7beb78", 
                                                  "Correct Model" = "#d55ae8",
                                                  "LASSO Model" = "#fc031c",
                                                  "Ridge Model" = "#ff983d")) +
          guides(fill = "none")
      }
    } else {
      if (num_legends >= 3) {
        plot <- plot + 
          scale_fill_manual(name = "", values = c("Selected Model" = "#5ab9e8", 
                                                  "Full Model" = "#7beb78", 
                                                  "Correct Model" = "#d55ae8"))
      } else {
        plot <- plot + 
          scale_fill_manual(name = "", values = c("Selected Model" = "#5ab9e8", 
                                                  "Full Model" = "#7beb78", 
                                                  "Correct Model" = "#d55ae8")) +
          guides(fill = "none")
      }
    }
    
    if (startsWith(beta, "Beta") && !is.null(Beta.seq)) {
      beta_number <- as.numeric(sub("Beta", "", beta))
      if (!is.na(beta_number) && beta_number + 1 <= length(Beta.seq)) {
        plot <- plot + geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed", color = "black")
      }
    }
    
    plot <- plot + 
      labs(title = plot_title, y = "Density", x = plot_x_label) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(size = 11),
            legend.position = "bottom")
    
    individual_plots[[length(individual_plots) + 1]] <- plot
  }
  
  # Add caption with the simulation settings
  caption_text <- paste("Simulation Settings:",
                        "Sample Size:", SampleSize, "|",
                        "Error Level:", Error, "|",
                        "Correlation Setting:", Correlation, "|",
                        "Information Criterion:", IC, "|",
                        "Model Selection Type:", MS.Type)
  
  combined_plot <- wrap_plots(individual_plots) + 
    plot_layout(ncol = 2, guides = "collect") &
    plot_annotation(title = "Coefficient Estimate Mixture Sampling Distributions",
                    caption = caption_text, 
                    theme = theme(plot.caption = element_text(size = 10, face = "italic", hjust = 0.5),
                                  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                  legend.position = c(0.9, 0.1),        # Positions legend inside bottom-right corner
                                  legend.justification = c(1, 0)        # Aligns legend's top-right corner to position
                                  )       
                    )        
  print(combined_plot)
  
  return(combined_plot)
}



# Confidence Interval Coverage Function
CI.Coverage.Fun <- function(Array, BetaTrue){
  # Empty object (add column for captured)
  CI.Array <- array(NA, dim = c(dim(Array)[1], dim(Array)[2] + 1, dim(Array)[3]))
  ModelNames <- rep(0, dim(Array)[3])
  
  for(i in 1:dim(Array)[3]){
    # Current Matrix
    CI.Mat <- Array[, , i]
    
    # Does the confidence interval encompass the true beta?
    Cap <- (BetaTrue >= CI.Mat[, 1]) & (BetaTrue <= CI.Mat[, 2])
    
    # Slap that bad boy into the array
    ArrayInput <- cbind(CI.Mat, Cap)
    CI.Array[, , i] <- ArrayInput 

    # Record the model 
    ModelNames[i] <- paste(rownames(CI.Mat)[!is.na(ArrayInput[, 3])][-1] , collapse = "")
  }
  
  # Name it
  dimnames(CI.Array) <- list(c("(Intercept)", paste("Beta", 1:(length(BetaTrue) - 1), sep = "")), 
                             c("Lower", "Upper", "Captured"),
                             ModelNames)
  
  # Calculate coverage across simulations
  p <- length(BetaTrue)
  Coverage.Mat <- matrix(data = NA, nrow = 2, ncol = p)
  
  for(i in 1:p){
    # Mean across ith row, 3rd column and all slices
    Coverage.Mat[1, i] <- mean(CI.Array[i, 3, ], na.rm = TRUE)
    
    # Count across ith row, 3rd column and all slices
    Coverage.Mat[2, i] <- sum(!is.na(CI.Array[i, 2, ]))
  }
  
  # Naming
  colnames(Coverage.Mat) <- c("(Intercept)", paste("Beta", 1:(length(BetaTrue) - 1), sep = ""))
  rownames(Coverage.Mat) <- c("Coverage", "Counts")
  
  # Output
  return(list(ConfidenceIntervals = CI.Array, Coverage.Mat = Coverage.Mat))
}

# Could put in model selection functions to have that in the output... was just applying it after the fact right now

# Post-Model Selection Coefficient vs Correct Model coefficients in a dataframe for tabulating
Model.Selection.Breakdown.Fun <- function(SimResults, MS.Type = c("Forward", "Backward", "Both", "LASSO", "Ridge"), 
                                          SampleSize = c("LowSS", "HighSS"), Error = c("LowError", "HighError"), 
                                          Correlation = c("Independent", "ModerateCorrelation"), IC = c("AIC", "BIC"),
                                          Model = c("Post.Model.Selection", "Correct.Model", "Full.Model"), 
                                          BetaSetting = c("PatternA", "PatternB"), Sim.Settings = TRUE) {
  ### Default Settings
  MS.Type <- match.arg(MS.Type)
  SampleSize <- match.arg(SampleSize)
  Error <- match.arg(Error)
  Correlation <- match.arg(Correlation)
  IC <- match.arg(IC)
  Model <- match.arg(Model)
  BetaSetting <- match.arg(BetaSetting)
  
  # Parameter Patterns
  if (BetaSetting == "PatternA") {
    True.Betas <- c(0, c(0.125, 0.25, 0.50, 1, 2))
  } else {
    True.Betas <- c(0, c(-1.5, -0.75, 0, 0.75, 1.5))
  }
  
  # Penalized or Sequential Model Selection
  if (MS.Type == "LASSO" || MS.Type == "Ridge") {
    TableData <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["LASSO.Ridge"]][[MS.Type]]
    Correct.Mod.Data <- NULL
    Full.Mod.Data <- NULL
    
  } else {
    TableData <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][[Model]]
    Correct.Mod.Data <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Correct.Model"]]
    Full.Mod.Data <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Full.Model"]]
  }
  
  # Load in all the fitted model coefficients
  AllMods <- SimResults[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["AllModels"]]
  UniRegs <- unique(AllMods[["Regressor"]])
  avg_all_mods <- data.frame(matrix(NA, nrow = length(UniRegs), ncol = dim(AllMods)[2]))
  for(i in 1:length(UniRegs)){
    TempData <- subset(AllMods, Regressor == UniRegs[i])
    avg_all_mods[i, ] <- c(TempData[1, 1], colMeans(TempData[,-1]))
  }
  avg_all_mods <- cbind.data.frame(avg_all_mods[,1], apply(avg_all_mods[,-1], 2, function(x) as.numeric(x)))
  colnames(avg_all_mods) <- colnames(AllMods) <- c("Regressors", paste("Beta", 0:(dim(AllMods)[2] - 2), "_E", sep = ""))
  
  # Tabulate Different Model Frequencies
  Reg.Table <- table(TableData[["Regressors"]])
  
  # Different Models Selected
  Mods <- names(Reg.Table)
  
  # Frequency of Each Model
  Mods.Freq <- as.vector(Reg.Table / nsim)
  
  # Subset the data for each model
  post.subsets <- list()
  for (mod in Mods) {
    # Subset based on model
    subs <- subset(TableData, Regressors == mod)
    
    # Remove rows that are not the coefficient values
    post.subsets[[mod]] <- subs[, startsWith(names(subs), "Beta")]
  }
  
  # Post-Model Selection Coefficients for each model
  PostCoefMeans <- t(sapply(post.subsets, colMeans, na.rm = TRUE))
  
  ### If LASSO / Ridge
  if(MS.Type == "LASSO" || MS.Type == "Ridge"){
    # Bias: True Beta - Post-Model Selection Coefficients
    if (is.null(True.Betas)) {
      stop("True.Betas must be provided to compute bias.")
    }
    
    # Make True.Betas a matrix with the same row names as PostCoefMeans
    True.Beta.Mat <- matrix(rep(True.Betas[-1], each = nrow(PostCoefMeans)), 
                            nrow = nrow(PostCoefMeans), byrow = FALSE)
    colnames(True.Beta.Mat) <- paste("Beta", 1:(length(True.Betas) - 1), sep = "")
    
    # Calculate bias
    Bias.Mat <- PostCoefMeans - True.Beta.Mat
    colnames(Bias.Mat) <- paste0(colnames(Bias.Mat), "_B")
    
    # Combine three matrices and order as desired
    Comb.Mat <- cbind(PostCoefMeans, Bias.Mat)
    
    # Extract base names (e.g., "Beta0" from "Beta0_C")
    base_names <- sub("(_.*)?$", "", colnames(Comb.Mat))
    
    # Order columns by numeric BetaX and then by suffix (_B, _C, _D, etc.)
    col_order <- order(as.numeric(sub("Beta", "", base_names)), colnames(Comb.Mat))
    
    # Reorder matrix
    Combined.Mat <- Comb.Mat[, col_order]
    
    # Get rid of row names since they'll be their own column
    rownames(Combined.Mat) <- NULL
    
    # Put back the model and model frequency columns
    if(is.null(dim(Combined.Mat))){
      Combined.df <- cbind.data.frame(Model = Mods, Chosen.Freq = Mods.Freq, t(Combined.Mat))
    } else {
      Combined.df <- cbind.data.frame(Model = Mods, Chosen.Freq = Mods.Freq, Combined.Mat)
    }
    
    # Add simulation settings to table as factor variables
    if (Sim.Settings == TRUE) {
      Combined.df["SampleSize"] <- SampleSize
      Combined.df["Error"] <- Error
      Combined.df["Correlation"] <- Correlation
      Combined.df["MS.Type"] <- MS.Type
      Combined.df["BetaSetting"] <- BetaSetting
    }

    # Output
    # cat("Simulation Settings:",
    #     "BetaSetting:" , BetaSetting, "|",
    #     "Sample Size:", SampleSize, "|",
    #     "Error Level:", Error, "|",
    #     "Correlation Setting:", Correlation, "|",
    #     "Information Criterion:", IC, "|",
    #     "Model Selection Type:", MS.Type, "\n")
    # print(Combined.df)
    return(Combined.df)
  }
  
  ### Correct model means
  # Append model list from post-model selection
  Cor.Mod.Data <- cbind.data.frame(Mods = TableData[["Regressors"]], Correct.Mod.Data)
  
  # Subset the data for each model
  correct.subsets <- list()
  for (mod in Mods) {
    # Subset based on model
    Subs <- subset(Cor.Mod.Data, Mods == mod)
    
    # Remove rows that are not the coefficient values
    correct.subsets[[mod]] <- Subs[, startsWith(names(Subs), "Beta")]
  }
  
  # Correct coefficient means
  CorrectCoefMeans <- t(sapply(correct.subsets, colMeans, na.rm = TRUE))
  
  # Change colnames so they don't match exactly
  colnames(CorrectCoefMeans) <- paste0(colnames(CorrectCoefMeans), "_C")
  
  # Find difference between two matrices
  DiffCoefMeans <- CorrectCoefMeans - PostCoefMeans
  
  # Change colnames again so they don't match exactly
  colnames(DiffCoefMeans) <- sub(".$", "D", colnames(CorrectCoefMeans))
  
  # Bias: True Beta - Post-Model Selection Coefficients
  if (is.null(True.Betas)) {
    stop("True.Betas must be provided to compute bias.")
  }

  # Make True.Betas a matrix with the same row names as PostCoefMeans
  True.Beta.Mat <- matrix(rep(True.Betas, each = nrow(PostCoefMeans)), 
                          nrow = nrow(PostCoefMeans), byrow = FALSE)
  colnames(True.Beta.Mat) <- colnames(PostCoefMeans)
  
  # Calculate bias
  Bias.Mat <- True.Beta.Mat - PostCoefMeans
  colnames(Bias.Mat) <- paste0(colnames(Bias.Mat), "_B")
  
  # Combine three matrices and order as desired
  Comb.Mat <- cbind(PostCoefMeans, CorrectCoefMeans, DiffCoefMeans, Bias.Mat)
  
  # Extract base names (e.g., "Beta0" from "Beta0_C")
  base_names <- sub("(_.*)?$", "", colnames(Comb.Mat))
  
  # Order columns by numeric BetaX and then by suffix (_B, _C, _D, etc.)
  col_order <- order(as.numeric(sub("Beta", "", base_names)), colnames(Comb.Mat))
  
  # Reorder matrix
  Combined.Mat <- Comb.Mat[, col_order]
  
  # Get rid of row names since they'll be their own column
  rownames(Combined.Mat) <- NULL
  
  # Put back the model and model frequency columns
  if(is.null(dim(Combined.Mat))){
    Combined.df <- cbind.data.frame(Model = Mods, Chosen.Freq = Mods.Freq, t(Combined.Mat))
  } else {
    Combined.df <- cbind.data.frame(Model = Mods, Chosen.Freq = Mods.Freq, Combined.Mat)
  }
  
  # Add in no model selection coefficient
  Combined.df <- merge(Combined.df, subset(avg_all_mods, Regressors %in% unique(Combined.df[["Model"]])),
                       by.x = "Model", by.y = "Regressors", all = TRUE)

  # Add simulation settings to table as factor variables
  if (Sim.Settings == TRUE) {
    Combined.df["SampleSize"] <- SampleSize
    Combined.df["Error"] <- Error
    Combined.df["Correlation"] <- Correlation
    Combined.df["IC"] <- IC
    Combined.df["MS.Type"] <- MS.Type
    Combined.df["BetaSetting"] <- BetaSetting
  }
  
  # Output
  # cat("Simulation Settings:",
  #     "BetaSetting:" , BetaSetting, "|",
  #     "Sample Size:", SampleSize, "|",
  #     "Error Level:", Error, "|",
  #     "Correlation Setting:", Correlation, "|",
  #     "Information Criterion:", IC, "|",
  #     "Model Selection Type:", MS.Type, "\n")
  # print(Combined.df)
  return(Combined.df)
}

# Reshaping smooshed tables into desired plotting dataframe
Reshape.Fun <- function(Smooshed.Table.Results){
  # Step 1: Pivot longer to stack all Beta and Bias columns
  df_long <- Smooshed.Table.Results %>%
    pivot_longer(
      cols = matches("^Beta\\d+(_B)?$"),  # Match Beta0, Beta1, ..., Beta0_B, etc.
      names_to = "BetaType",
      values_to = "Value"
    )
  
  # Step 2: Split into 'Beta' and 'Type' (Estimate or Bias)
  df_long <- df_long %>%
    mutate(
      Type = if_else(str_detect(BetaType, "_B$"), "Bias", "Estimate"),
      Beta = str_remove(BetaType, "_B$")
    ) %>%
    select(-BetaType)
  
  # Step 3: Add true beta values 
  df_longA <- subset(df_long, BetaSetting == "PatternA")
  True.BetasA <- c(0, 0.125, 0.25, 0.50, 1, 2)
  
  true_betasA <- tibble(
    Beta = paste0("Beta", 0:(length(True.BetasA) - 1)),
    True = True.BetasA
  )
  
  df_longA <- df_longA %>%
    left_join(true_betasA, by = "Beta") %>%
    relocate(True, .after = Value)
  
  df_longB <- subset(df_long, BetaSetting == "PatternB")
  True.BetasB <- c(0, -1.5, -0.75, 0, 0.75, 1.5)
  
  true_betasB <- tibble(
    Beta = paste0("Beta", 0:(length(True.BetasB) - 1)),
    True = True.BetasB
  )
  
  df_longB <- df_longB %>%
    left_join(true_betasB, by = "Beta") %>%
    relocate(True, .after = Value)
  
  df_long <- bind_rows(df_longA, df_longB)
  
  
  # Step 4: Pivot wider so Estimate and Bias are their own columns
  df_long_wider <- df_long[2 * c(0:(nrow(df_long) / 2)) + 1, ]
  df_long_wider[["Estimate"]] <- df_long_wider$Value
  df_long_wider[["Bias"]] <- df_long[2 * c(0:(nrow(df_long) / 2)) + 2, ][["Value"]]
  df_long_wider <- subset(df_long_wider, select = (-c(Type, Value)))
  PlotData <- as.data.frame(df_long_wider)
  
  return(PlotData)
}
