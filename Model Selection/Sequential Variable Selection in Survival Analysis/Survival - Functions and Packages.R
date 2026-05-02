### Functions and Packages
## Load packages
library(survival)
library(tidyverse)

# Data generation function
Surv.Data.Gen.Fun <- function(n, Beta, Sigma, Treat = TRUE, prob = 0.5,
                              CensDist = "exponential",
                              CensPar = list(rate = 0.5),
                              ObsDist = "exponential",
                              ObsPar = list(rate = 0.5, shape = 1), ...){
  # Number of predictors
  p <- length(Beta) - 1
  
  if(Treat == TRUE) {
    # Generate predictor matrix X with specified covariance structure Sigma along with a treatment variable
    X <- cbind(rbinom(n, size = 1, prob = prob),
               MASS::mvrnorm(n = n, mu = rep(0, p-1), Sigma = Sigma))
    
  } else {
    # Generate predictor matrix X with specified covariance structure Sigma
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  }
  
  # Observed times based on the specified distribution
  if (ObsDist == "exponential") {
    # Survreg negates the signs in the fitted models
    rate <- ObsPar[["rate"]] * exp(-Beta[1] + X %*% -Beta[-1])
    Observations <- rexp(n, rate = rate)
    
  } else if(ObsDist == "lognormal"){
    meanlog <- Beta[1] + X %*% Beta[-1]
    sdlog <- ObsPar[["shape"]]
    Observations <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
    
  } else if(ObsDist == "weibull"){
    shape <- ObsPar[["shape"]]
    scale <- ObsPar[["scale"]] * exp(Beta[1] + X %*% Beta[-1])
    Observations <- rweibull(n, shape = shape, scale = scale)
    
  } else if(ObsDist == "Cox" | ObsDist == "cox"){
    # Based Gompertz distribution
    shape <- ObsPar[["shape"]]
    scale <- ObsPar[["scale"]]
    Observations <- log(1 - log(runif(n))/(exp(Beta[1] + X %*% Beta[-1])*shape))/scale
    
  } else {
    stop("Other observation distributions not available")
    
  }
  
  # Censoring distribution
  if(CensDist == "exponential"){
    Censored.Vec <- rexp(n, CensPar[["rate"]])
    
  } else {
    stop("Other censoring distributions not available")
    
  }
  
  # Observed times
  t <- pmin(Observations, Censored.Vec)
  
  # Event indicator
  delta <- as.numeric(Observations <= Censored.Vec)
  
  return(data.frame(time = t, delta, X))
}

# Sequential model selection function for survival models
Survival.MS.Fun <- function(Data, Beta, Model = "Cox", MS.Type = "forward", IC = "AIC",
                            fullscope = as.formula("Surv(time, delta) ~ ."),
                            nullscope = as.formula("Surv(time, delta) ~ X1"), # Attempting to force in treatment
                            dist = "weibull") {
  # Housekeeping
  da.data <- as.data.frame(Data)
  colnames(da.data) <- c("time", "delta", paste("X", 1:(ncol(da.data) - 2), sep = ""))
  
  # Perform Sequential Selection based on AIC
  if (IC == "AIC") {
    if(Model == "Cox"){
      # Null Cox model (only baseline hazard)
      Null.Cox <- coxph(as.formula(nullscope), data = da.data)
      
      # Full Cox model (all variables)
      Full.Cox <- coxph(as.formula(fullscope), data = da.data)
      
      # Forward Selection
      if(MS.Type == "forward"){
        # Cox sequential selection
        Sel.Mod <- step(Null.Cox, scope = list(lower = Null.Cox, upper = Full.Cox), 
                        direction = MS.Type, k = 2, trace = 0)
        
        # Backwards or Bidirectional Selection
      } else if (MS.Type == "backward" | MS.Type == "both") {
        # Cox sequential selection
        Sel.Mod <- step(Full.Cox, scope = list(lower = Null.Cox, upper = Full.Cox), 
                        direction = MS.Type, k = 2, trace = 0)
        
      } else {
        stop("Model selection method not supported")
        
      }
      
    } else if (Model == "AFT") {
      # Null AFT model (only intercept)
      Null.AFT <- survreg(as.formula(nullscope), data = da.data, dist = dist, beta = Beta,
                          control = survreg.control(maxiter = 500))
      
      # Full AFT model (all variables)
      Full.AFT <- survreg(as.formula(fullscope), data = da.data, dist = dist, beta = Beta, 
                          control = survreg.control(maxiter = 500))
      
      # Forward Selection
      if(MS.Type == "forward"){
        # AFT sequential selection
        Sel.Mod <- step(Null.AFT, scope = list(lower = Null.AFT, upper = Full.AFT), 
                        direction = MS.Type, k = 2, trace = 0)
        
        # Backwards or Bidirectional Selection
      } else if (MS.Type == "backward" | MS.Type == "both") {
        # AFT sequential selection
        Sel.Mod <- step(Full.AFT, scope = list(lower = Null.AFT, upper = Full.AFT),
                        direction = MS.Type, k = 2, trace = 0)

      } else {
        stop("Model selection method not supported")
        
      }
    }
    # Perform Sequential Selection based on BIC  
  } else if (IC == "BIC") {
    if(Model == "Cox"){
      # Null Cox model (only baseline hazard)
      Null.Cox <- coxph(as.formula(nullscope), data = da.data)
      
      # Full Cox model (all variables)
      Full.Cox <- coxph(as.formula(fullscope), data = da.data)
      
      # Forward Selection
      if(MS.Type == "forward"){
        # Cox sequential selection
        Sel.Mod <- step(Null.Cox, scope = list(lower = Null.Cox, upper = Full.Cox), 
                        direction = MS.Type, k = log(nrow(da.data)), trace = 0)
        
        # Backwards or Bidirectional Selection
      } else if (MS.Type == "backward" | MS.Type == "both") {
        # Cox sequential selection
        Sel.Mod <- step(Full.Cox, scope = list(lower = Null.Cox, upper = Full.Cox),
                        direction = MS.Type, k = log(nrow(da.data)), trace = 0)
        
      } else {
        stop("Model selection method not supported")
        
      }
      
    } else if (Model == "AFT") {
      # Null AFT model (only intercept)
      Null.AFT <- survreg(as.formula(nullscope), data = da.data, dist = dist, beta = Beta, 
                          control = survreg.control(maxiter = 500))
      
      # Full AFT model (all variables)
      Full.AFT <- survreg(as.formula(fullscope), data = da.data, dist = dist, beta = Beta, 
                          control = survreg.control(maxiter = 500))
      
      # Forward Selection
      if(MS.Type == "forward"){
        # AFT sequential selection
        Sel.Mod <- step(Null.AFT, scope = list(lower = Null.AFT, upper = Full.AFT), 
                        direction = MS.Type, k = log(nrow(da.data)), trace = 0)
        
        # Backwards or Bidirectional Selection
      } else if (MS.Type == "backward" | MS.Type == "both") {
        # AFT sequential selection
        Sel.Mod <- step(Full.AFT, scope = list(lower = Null.AFT, upper = Full.AFT),
                        direction = MS.Type, k = log(nrow(da.data)), trace = 0)
        
      } else {
        stop("Model selection method not supported")
        
      }
    }
    
  } else {
    stop("Information criteria not supported")
  }
  
  return(Sel.Mod)
}

# Function to apply column names based on matrix dimensions
ApplyColnames <- function(mat) {
  if (is.matrix(mat)) {
    n_cols <- ncol(mat)
    
    if (n_cols == length(Beta)) {
      colnames(mat) <- c("Intercept", paste("Beta", 1:(n_cols - 1), sep = ""))
      
    } else if (n_cols == 4) {
      colnames(mat) <- c("AFT1", "AFT2", "AFT3", "Cox")
      
    } else if (n_cols == 2) {
      colnames(mat) <- c("Concord", "Concord SE")
    }
  }
  return(mat)
}


### Helper Function for the function below 
# Does both different types of averaging preserving the regressor and distribution combinations
# Also notes the count of models 
average_by_group <- function(df, na_as_zero = FALSE) {
  
  # Get unique combinations
  regressors_list <- unique(df$Regressors)
  distribution_list <- unique(df$Distribution)
  
  result_list <- list()
  
  for (r in regressors_list) {
    for (d in distribution_list) {
      # Subset the data
      temp <- df[df$Regressors == r & df$Distribution == d, ]
      
      if (nrow(temp) > 0) {
        # Select only numeric columns
        numeric_cols <- sapply(temp, is.numeric)
        temp_numeric <- temp[, numeric_cols, drop = FALSE]
        
        # Handle NAs
        # Treat NAs as zero
        if (na_as_zero) {
          temp_numeric[is.na(temp_numeric)] <- 0
          temp_means <- colMeans(temp_numeric, na.rm = FALSE)
          # Ignore NAs in averaging
        } else {
          temp_means <- colMeans(temp_numeric, na.rm = TRUE)
        }
        
        # Rebuild the row
        temp_result <- data.frame(as.list(temp_means))
        
        # Add Frequency / count column (number of rows being averaged)
        temp_result$Count <- nrow(temp)
        
        # Add Regressors and Distribution
        temp_result$Regressors <- r
        temp_result$Distribution <- d
        
        # Store result
        result_list[[length(result_list) + 1]] <- temp_result
      }
    }
  }
  
  # Combine into one data frame
  final_result <- do.call(rbind, result_list)
  
  # Clean up rownames
  rownames(final_result) <- NULL
  
  # Return
  return(final_result)
}

# Adds captured column to arrays and calculates average coverage
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

# Helper function for the agreement table
Pull_Regs <- function(df, Scenario) {
  if (Scenario == "Scenario11") {
    df1 <- df[["AFT1_Weibull"]][["MS_Coef"]][["Regressors"]]
    df2 <- df[["AFT2_Exponential"]][["MS_Coef"]][["Regressors"]]
    df3 <- df[["AFT3_Lognormal"]][["MS_Coef"]][["Regressors"]]
  } else if (Scenario == "Scenario12") {
    df1 <- df[["AFT1_Weibull"]][["MS_Sel_Coef"]][["Regressors"]]
    df2 <- df[["AFT2_Exponential"]][["MS_Sel_Coef"]][["Regressors"]]
    df3 <- df[["AFT3_Lognormal"]][["MS_Sel_Coef"]][["Regressors"]]
  } else if (Scenario == "Scenario21") {
    df1 <- df[["AFT1_Weibull"]][["Sel_Coef"]][["Regressors"]]
    df2 <- df[["AFT2_Exponential"]][["Sel_Coef"]][["Regressors"]]
    df3 <- df[["AFT3_Lognormal"]][["Sel_Coef"]][["Regressors"]]
  } else if (Scenario == "Scenario22") {
    df1 <- df[["AFT1_Weibull"]][["Sel_MS_Coef"]][["Regressors"]]
    df2 <- df[["AFT2_Exponential"]][["Sel_MS_Coef"]][["Regressors"]]
    df3 <- df[["AFT3_Lognormal"]][["Sel_MS_Coef"]][["Regressors"]]
  }
  
  # Combine while checking for character "NA" and empty strings
  Regs <- ifelse(df1 != "NA" & df1 != "", df1, 
                 ifelse(df2 != "NA" & df2 != "", df2, df3))
  
  return(Regs)
}

# Function that tabulates simulations for each unique simulation setting combination
SimulationTabulation.Fun <- function(SimResults, CI, correlation = c("Independent", "ModerateCorrelation"),
                                     censoring = c("LowCens", "HighCens"),
                                     distribution = c("weibull", "exponential", "lognormal", "cox"),
                                     BetaTruth = Beta) {
  ### Default Settings
  correlation <- match.arg(correlation)
  censoring <- match.arg(censoring)
  distribution <- match.arg(distribution)
  
  ### Pulling desired output based on simulation setting
  TableData <- SimResults[[correlation]][[censoring]][[distribution]]
  CIData <- CI[[correlation]][[censoring]][[distribution]]
  
  ### Just keeping list elements with coefficient values
  CoefTable <- list(AFT1_Weibull = TableData[["AFT1_Weibull"]][grepl("_Coef$", names(TableData[["AFT1_Weibull"]]))],
                    AFT2_Exponential = TableData[["AFT2_Exponential"]][grepl("_Coef$", names(TableData[["AFT2_Exponential"]]))],
                    AFT3_Lognormal = TableData[["AFT3_Lognormal"]][grepl("_Coef$", names(TableData[["AFT3_Lognormal"]]))],
                    Cox_Regression = TableData[["Cox_Regression"]][grepl("_Coef$", names(TableData[["Cox_Regression"]]))])

  ### Creating a regressor column for each element
  CoefTable <- lapply(CoefTable, FUN = function(liss){
    liss <- lapply(liss, FUN = function(df){
      df <- as.data.frame(df)
      df[["Regressors"]] <- rep(NA, nrow(df))
      for(i in 1:nrow(df)){
        # If only the first element is present, label it the intercept model
        if (!is.na(df[i, ][1]) && all(is.na(df[i, ][-1]))) {
          df[i, "Regressors"] <- "Int"
        
        # If all NA, then it's an NA
        } else if (sum(!is.na(df[i, ])) == 0) {
          df[i, "Regressors"] <- "NA"
        
        # For all other cases give a letter to each present beta
        } else {
          df[i, "Regressors"] <- paste(LETTERS[1:(ncol(df) - 2)][!is.na(df[i, -1])], collapse = "")
        }
      }
      return(df)
      })
  })
  
  # Aggregation Table
  ChosenModelDF <- SimResults[[correlation]][[censoring]][[distribution]][["Miscellaneous"]][["Chosen.Model"]]
  ChosenModelDF <- data.frame(sapply(rep(ChosenModelDF, each = 2), function(x) as.vector(x)))
  colnames(ChosenModelDF) <- paste("Model.", c("Scenario11", "Scenario12", "Scenario21", "Scenario22"), sep = "")
  
  # Pulling Regressors
  ChosenModelDF <- data.frame(append(ChosenModelDF, list(Regs.Scenario11 = Pull_Regs(CoefTable, Scenario = "Scenario11")), after = 1))
  ChosenModelDF <- data.frame(append(ChosenModelDF, list(Regs.Scenario12 = Pull_Regs(CoefTable, Scenario = "Scenario12")), after = 3))
  ChosenModelDF <- data.frame(append(ChosenModelDF, list(Regs.Scenario21 = Pull_Regs(CoefTable, Scenario = "Scenario21")), after = 5))
  ChosenModelDF <- data.frame(append(ChosenModelDF, list(Regs.Scenario22 = Pull_Regs(CoefTable, Scenario = "Scenario22")), after = 7))

  # Flatten lists
  CoefList <- do.call(c, CoefTable)
  CIList <- do.call(c, CIData)
  
  # Convert CI arrays into data frames
  CIList_df <- lapply(CIList, function(arr) {
    # Apply CI Coverage function - Really just to add the captured column (I'm lazy)
    Arr <- CI.Coverage.Fun(arr, BetaTrue = BetaTruth)[[1]]
    
    # Column names: Intercept_Lower, Intercept_Upper, Intercept_Captured ..., Beta5_Captured
    ci_colnames <- as.vector(t(outer(c("Intercept", paste0("Beta", 1:(length(BetaTruth) - 1))),
                                     c("Lower", "Upper", "Captured"), paste, sep = "_")))
    # Shape array into dataframe
    n <- dim(Arr)[3]
    slice_list <- lapply(1:n, function(i) {
      mat <- Arr[,,i]         # 6x3 matrix
      row <- as.vector(t(mat)) # Flatten row-wise
      setNames(as.data.frame(t(row)), ci_colnames)
    })
    
    return(do.call(rbind, slice_list))
  })
  
  # Final merge: cbind each df
  Merged <- mapply(cbind, CoefList, CIList_df, SIMPLIFY = FALSE)
  
  # Breaking up depending on which stage of the design we're at
  Full.ModelList <- Merged[grepl("\\.Full_Coef$", names(Merged))]
  MS.Merged <- Merged[grepl("\\.MS_Coef$", names(Merged))]
  MS.Sel.Merged <- Merged[grepl("\\.MS_Sel_Coef$", names(Merged))]
  Sel.Merged <- Merged[grepl("\\.Sel_Coef$", names(Merged))]
  Sel.MS.Merged <- Merged[grepl("\\.Sel_MS_Coef$", names(Merged))]
  
  # Collapsing each into one their own big data frame and getting rid of any rows that are all NAs
  Full.Modeldf <- bind_rows(Full.ModelList)
  Full.Modeldf[["Distribution"]] <- rep(sub("\\..*$", "", names(Full.ModelList)), each = nrow(Full.ModelList[[1]]))
  
  MS.Coefdf <- bind_rows(MS.Merged)
  MS.Coefdf[["Distribution"]] <- rep(sub("\\..*$", "", names(MS.Merged)), each = nrow(MS.Merged[[1]]))
  MS.Coefdf <- MS.Coefdf[!apply(subset(MS.Coefdf, select = -c(Regressors, Distribution)), 1, function(x) all(is.na(x))), ]
  
  MS.Sel.Coefdf <- bind_rows(MS.Sel.Merged)
  MS.Sel.Coefdf[["Distribution"]] <- rep(sub("\\..*$", "", names(MS.Sel.Merged)), each = nrow(MS.Sel.Merged[[1]]))
  MS.Sel.Coefdf <- MS.Sel.Coefdf[!apply(subset(MS.Sel.Coefdf, select = -c(Regressors, Distribution)), 1, function(x) all(is.na(x))), ]
  
  Sel.Coefdf <- bind_rows(Sel.Merged)
  Sel.Coefdf[["Distribution"]] <- rep(sub("\\..*$", "", names(Sel.Merged)), each = nrow(Sel.Merged[[1]]))
  Sel.Coefdf <- Sel.Coefdf[!apply(subset(Sel.Coefdf, select = -c(Regressors, Distribution)), 1, function(x) all(is.na(x))), ]
  
  Sel.MS.Coefdf <- bind_rows(Sel.MS.Merged)
  Sel.MS.Coefdf[["Distribution"]] <- rep(sub("\\..*$", "", names(Sel.MS.Merged)), each = nrow(Sel.MS.Merged[[1]]))
  Sel.MS.Coefdf <- Sel.MS.Coefdf[!apply(subset(Sel.MS.Coefdf, select = -c(Regressors, Distribution)), 1, function(x) all(is.na(x))), ]
  
  
  ### Conditional and Unconditional Averaging (Using helper function defined before)
  # Full model results (No NAs to worry about)
  Avg.Full.Modeldf <- data.frame(average_by_group(Full.Modeldf, na_as_zero = FALSE),
                                       Averaging = "NoNAs", Scenario = "S00")
  
  # Scenario 1 - Part 1 (No NAs to worry about)
  Avg.MS.Coefdf <- data.frame(average_by_group(MS.Coefdf, na_as_zero = FALSE),
                                    Averaging = "NoNAs", Scenario = "S11")
  
  # Scenario 1 - Part 2
  Avg.Cond.MS.Sel.Coefdf <- data.frame(average_by_group(MS.Sel.Coefdf, na_as_zero = FALSE),
                                             Averaging = "Conditional", Scenario = "S12") # Conditional
  Avg.Uncon.MS.Sel.Coefdf <- data.frame(average_by_group(MS.Sel.Coefdf, na_as_zero = TRUE),
                                              Averaging = "Unconditional", Scenario = "S12") # Unconditional
  
  # Scenario 2 - Part 1
  Avg.Cond.Sel.Coefdf <- data.frame(average_by_group(Sel.Coefdf, na_as_zero = FALSE),
                                          Averaging = "Conditional", Scenario = "S21") # Conditional
  Avg.Uncon.Sel.Coefdf <- data.frame(average_by_group(Sel.Coefdf, na_as_zero = TRUE),
                                           Averaging = "Unconditional", Scenario = "S21") # Unconditional
  
  # Scenario 2 - Part 2
  Avg.Cond.Sel.MS.Coefdf <- data.frame(average_by_group(Sel.MS.Coefdf, na_as_zero = FALSE),
                                             Averaging = "Conditional", Scenario = "S22") # Conditional
  Avg.Uncon.Sel.MS.Coefdf <- data.frame(average_by_group(Sel.MS.Coefdf, na_as_zero = TRUE), 
                                        Averaging = "Unconditional", Scenario = "S22") # Unconditional
  
  ### Slapping all into a dataframe and adding the simulation settings to said dataframe
  Final.df <- bind_rows(Avg.Full.Modeldf, 
                         Avg.MS.Coefdf,
                         Avg.Cond.MS.Sel.Coefdf,
                         Avg.Uncon.MS.Sel.Coefdf,
                         Avg.Cond.Sel.Coefdf,
                         Avg.Uncon.Sel.Coefdf,
                         Avg.Cond.Sel.MS.Coefdf,
                         Avg.Uncon.Sel.MS.Coefdf)
  
  # Adding current simulation settings
  Final.df["Correlation"] <- correlation
  Final.df["Censoring"] <- censoring
  Final.df["Gen.Distribution"] <- distribution
  
  # Output
  return(list(Final.df, ChosenModelDF))
}

# To reshape for plotting...
Reshape.Fun <- function(df) {
  # Step 0: Add row ID to track original rows
  df <- df %>% mutate(RowID = row_number())
  
  # Step 1: Rename Intercept to Beta0
  df <- df %>% rename_with(~ str_replace_all(., "^Intercept", "Beta0"))
  
  beta_names <- paste0("Beta", 1:4)
  
  # Step 2: Pivot estimates, CI, captured
  estimate_data <- df %>%
    select(RowID, all_of(beta_names)) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "Estimate")
  
  lower_data <- df %>%
    select(RowID, all_of(paste0(beta_names, "_Lower"))) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "Lower") %>%
    mutate(Beta = str_remove(Beta, "_Lower"))
  
  upper_data <- df %>%
    select(RowID, all_of(paste0(beta_names, "_Upper"))) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "Upper") %>%
    mutate(Beta = str_remove(Beta, "_Upper"))
  
  captured_data <- df %>%
    select(RowID, all_of(paste0(beta_names, "_Captured"))) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "Captured") %>%
    mutate(Beta = str_remove(Beta, "_Captured"))
  
  # Step 3: Pivot True and Bias
  true_data <- df %>%
    select(RowID, matches("^Beta_[1-4]_True$")) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "True") %>%
    mutate(Beta = paste0("Beta", str_extract(Beta, "[1-4]")))
  
  bias_data <- df %>%
    select(RowID, matches("^Beta_[1-4]_Bias$")) %>%
    pivot_longer(-RowID, names_to = "Beta", values_to = "Bias") %>%
    mutate(Beta = paste0("Beta", str_extract(Beta, "[1-4]")))
  
  # Step 4: Join all beta-level pieces by RowID + Beta
  beta_info <- estimate_data %>%
    left_join(lower_data, by = c("RowID", "Beta")) %>%
    left_join(upper_data, by = c("RowID", "Beta")) %>%
    left_join(captured_data, by = c("RowID", "Beta")) %>%
    left_join(true_data, by = c("RowID", "Beta")) %>%
    left_join(bias_data, by = c("RowID", "Beta"))
  
  # Step 5: Add back metadata
  meta_cols <- df %>%
    select(-matches("^(Beta[1-4](_(Lower|Upper|Captured))?$|Beta_[1-4]_(True|Bias))"), -RowID)
  
  final_df <- beta_info %>%
    left_join(df %>% select(RowID), by = "RowID") %>%
    bind_cols(meta_cols[rep(1:nrow(meta_cols), each = 4), ]) %>%
    select(-RowID)
  
  return(final_df)
}
# Not general - only works for beta 1 to 4

### Confusion Matrix Function
# Can create matrix of just Model dist agreement and / or model regressor agreement
# Can match with chosen distribution if doing both variables (not working as intended / useless) - otherwise creates a composite variable
# Input NULL if not interested in variable
ConfusedTable.Fun <- function(df_list, ScenarioA, ScenarioB, var1 = "Model", var2 = "Reg", match_dist = FALSE) {
  matrix_list <- list()
  for(i in 1:length(df_list)){
    # Split into desired data frames
    df1 <- select(df_list[[i]], matches(ScenarioA))
    df2 <- select(df_list[[i]], matches(ScenarioB))
    
    if(is.null(var1)){
      var2_df1 <- select(df1, matches(var2))
      var2_df2 <- select(df2, matches(var2))
      
      # Equalize lengths by padding with NA
      length(var2_df1) <- length(var2_df2) <- max(length(var2_df1), length(var2_df2))
      
      # Create the confusion matrix
      confusion_matrix <- table(unlist(var2_df1), unlist(var2_df2), useNA = "ifany")
      
      # Confused list
      matrix_list[[i]] <- confusion_matrix
      
      next
    }
    
    # Extract the model chosen
    var1_df1 <- select(df1, matches(var1))
    var1_df2 <- select(df2, matches(var1))
    
    # Check if the second variable is provided and exists
    if (!is.null(var2)) {
      var2_df1 <- select(df1, matches(var2))
      var2_df2 <- select(df2, matches(var2))
      
      # Handle matching distributions if requested
      if (match_dist) {
        matching_rows <- var1_df1 == var1_df2 & !is.na(var1_df1) & !is.na(var1_df2)
        var1_df1 <- var2_df1[matching_rows]
        var1_df2 <- var2_df2[matching_rows]
      } else {
        # Create composite variables (e.g., "Regressor_Distribution")
        var1_df1 <- paste(unlist(var1_df1), unlist(var2_df1), sep = "_")
        var1_df2 <- paste(unlist(var1_df2), unlist(var2_df2), sep = "_")
      }
    }
    
    if(is.null(var2)){
      var1_df1 <- unlist(var1_df1)
      var1_df2 <- unlist(var1_df2)
    }
    
    # Equalize lengths by padding with NA
    length(var1_df1) <- length(var1_df2) <- max(length(var1_df1), length(var1_df2))
    
    # Create the confusion matrix
    confusion_matrix <- table(var1_df1, var1_df2, useNA = "ifany")
    
    # Confused list
    matrix_list[[i]] <- confusion_matrix
  }
  
  # Step 1: Find all unique row and column names
  all_rows <- unique(unlist(lapply(matrix_list, rownames)))
  all_cols <- unique(unlist(lapply(matrix_list, colnames)))
  
  # Step 2: Initialize a result matrix filled with zeros
  result_matrix <- matrix(0, nrow = length(all_rows), ncol = length(all_cols), 
                          dimnames = list(all_rows, all_cols))
  
  # Step 3: Sum matrices into the result matrix
  for (mat in matrix_list) {
    # Match rows and columns to the result matrix positions
    row_idx <- match(rownames(mat), all_rows)
    col_idx <- match(colnames(mat), all_cols)
    
    # Add the current matrix values to the result matrix
    result_matrix[row_idx, col_idx] <- result_matrix[row_idx, col_idx] + as.matrix(mat, na.rm = TRUE)
  }
  
  return(result_matrix)
}

### Testing
# ConfusedTable.Fun(Chosen.Mod.Results, ScenarioA = "Scenario11", ScenarioB = "Scenario12")

##############################
# Rewrite of survreg to allow flexible initialization

survreg <- function (formula, data, weights, subset, na.action, dist = "weibull",
          init = NULL, beta=NULL, scale = 0, control, parms = NULL, model = FALSE,
          x = FALSE, y = TRUE, robust = FALSE, cluster, score = FALSE,
          ...)
{
  Call <- match.call()
  if (missing(formula))
    stop("a formula argument is required")
  ss <- c("cluster", "offset")
  Terms <- if (missing(data))
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)

  if(is.null(beta) & is.null(init)){
    init = rep(0, length(attr(terms(formula, data=data),
                              "term.labels"))+1)
  }else if(is.null(init)){
    init = 0.8*c(beta[1], beta[which(names(data)[-(1:2)] %in% (attr(terms(formula, data=data),
                            "term.labels")))+1])
  }
  
  #print(names(data))
  #print(attr(terms(formula, data=data), "term.labels"))
  #print(init)

  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1)
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1))
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0)
      stop("formula must have a Surv response")
    if (is.null(Call$cluster))
      Call$cluster <- attr(Terms, "variables")[[1 + tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    Terms <- drop.special(Terms, tcl)
    formula <- Call$formula <- formula(Terms)
  }
  indx <- match(c("formula", "data", "weights", "subset", "na.action",
                  "cluster"), names(Call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  special <- c("strata")
  temp$formula <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  m <- eval(temp, parent.frame())
  Terms <- attr(m, "terms")
  weights <- model.extract(m, "weights")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type == "counting")
    stop("start-stop type Surv objects are not supported")
  if (type == "mright" || type == "mcounting")
    stop("multi-state survival is not supported")
  cluster <- model.extract(m, "cluster")
  if (length(cluster)) {
    if (missing(robust))
      robust <- TRUE
    cluster <- as.numeric(as.factor(cluster))
  }
  else if (robust)
    cluster <- 1:nrow(Y)
  strats <- attr(Terms, "specials")$strata
  dropx <- NULL
  if (length(strats)) {
    temp <- untangle.specials(Terms, "strata", 1)
    dropx <- temp$terms
    if (length(temp$vars) == 1)
      strata.keep <- m[[temp$vars]]
    else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
    strata <- as.numeric(strata.keep)
    nstrata <- max(strata)
  }
  else {
    nstrata <- 1
    strata <- 0
  }
  if (length(dropx)) {
    newTerms <- Terms[-dropx]
    attr(newTerms, "intercept") <- attr(Terms, "intercept")
  }
  else newTerms <- Terms
  X <- model.matrix(newTerms, m)
  assign <- lapply(attrassign(X, newTerms)[-1], function(x) x -
                     1)
  xlevels <- .getXlevels(newTerms, m)
  contr.save <- attr(X, "contrasts")
  if (!all(is.finite(X)))
    stop("data contains an infinite predictor")
  n <- nrow(X)
  nvar <- ncol(X)
  offset <- model.offset(m)
  if (length(offset) == 0 || all(offset == 0))
    offset <- rep(0, n)
  if (is.character(dist)) {
    dist <- match.arg(dist, names(survreg.distributions))
    dlist <- survreg.distributions[[dist]]
    if (is.null(dlist))
      stop(paste(dist, ": distribution not found"))
  }
  else if (is.list(dist))
    dlist <- dist
  else stop("Invalid distribution object")
  if (!survregDtest(dlist))
    stop("Invalid distribution object")
  logcorrect <- 0
  Ysave <- Y
  if (!is.null(dlist$trans)) {
    tranfun <- dlist$trans
    exactsurv <- Y[, ncol(Y)] == 1
    if (any(exactsurv)) {
      if (is.null(weights))
        logcorrect <- sum(log(dlist$dtrans(Y[exactsurv,
                                             1])))
      else logcorrect <- sum(weights[exactsurv] * log(dlist$dtrans(Y[exactsurv,
                                                                     1])))
    }
    if (type == "interval") {
      if (any(Y[, 3] == 3))
        Y <- cbind(tranfun(Y[, 1:2]), Y[, 3])
      else Y <- cbind(tranfun(Y[, 1]), Y[, 3])
    }
    else if (type == "left")
      Y <- cbind(tranfun(Y[, 1]), 2 - Y[, 2])
    else Y <- cbind(tranfun(Y[, 1]), Y[, 2])
    if (!all(is.finite(Y)))
      stop("Invalid survival times for this distribution")
  }
  else {
    if (type == "left")
      Y[, 2] <- 2 - Y[, 2]
    else if (type == "interval" && all(Y[, 3] < 3))
      Y <- Y[, c(1, 3)]
  }
  if (!is.null(dlist$scale)) {
    if (!missing(scale))
      warning(paste(dlist$name, "has a fixed scale, user specified value ignored"))
    scale <- dlist$scale
  }
  if (!is.null(dlist$dist))
    if (is.atomic(dlist$dist))
      dlist <- survreg.distributions[[dlist$dist]]
  else dlist <- dlist$dist
  ptemp <- dlist$parms
  if (is.null(ptemp)) {
    if (!is.null(parms))
      stop(paste(dlist$name, "distribution has no optional parameters"))
  }
  else {
    if (!is.numeric(ptemp))
      stop("Default parameters must be a numeric vector")
    if (!missing(parms)) {
      temp <- unlist(parms)
      indx <- match(names(temp), names(ptemp))
      if (any(is.na(indx)))
        stop("Invalid parameter names")
      ptemp[names(ptemp)] <- temp
    }
    parms <- ptemp
  }
  if (missing(control))
    control <- survreg.control(...)
  else control <- do.call("survreg.control", control)
  if (any(scale < 0))
    stop("Invalid scale value")
  if (any(scale > 0) && nstrata > 1)
    stop("The scale argument is not valid with multiple strata")
  pterms <- sapply(m, inherits, "coxph.penalty")
  if (any(pterms)) {
    if (any(grepl("frailty", names(pterms))))
      stop("survreg does not support frailty terms")
    pattr <- lapply(m[pterms], attributes)
    temp <- c(attr(Terms, "response"), attr(Terms, "offset"))
    if (length(dropx))
      temp <- c(temp, dropx + 1)
    pterms <- pterms[-temp]
    temp <- match((names(pterms))[pterms], attr(Terms, "term.labels"))
    ord <- attr(Terms, "order")[temp]
    if (any(ord > 1))
      stop("Penalty terms cannot be in an interaction")
    assign <- attrassign(X, newTerms)
    pcols <- assign[match(names(pterms[pterms]), names(assign))]
    fit <- survpenal.fit(X, Y, weights, offset, init = init,
                         controlvals = control, dist = dlist, scale = scale,
                         strata = strata, nstrat = nstrata, pcols, pattr,
                         parms = parms, assign)
  }
  else{
    n.init <- 1
    while(n.init < 10){
      fit <- survreg.fit(X, Y, weights, offset, init = init,
                         controlvals = control, dist = dlist, scale = scale, nstrat = nstrata,
                         strata, parms = parms)
      if (scale == 0) {
          nvar <- length(fit$coefficients) - nstrata
          fit$scale <- exp(fit$coefficients[-(1:nvar)])
          if (nstrata == 1)
            names(fit$scale) <- NULL
          else names(fit$scale) <- levels(strata.keep)
          fit$coefficients <- fit$coefficients[1:nvar]
          fit$idf <- 1 + nstrata
      }
      else {
        fit$scale <- scale
        fit$idf <- 1
      }
        
      singular <- (diag(fit$var) == 0)[1:length(fit$coefficients)]
      if (!any(singular)) {
        fit <- survreg.fit(X, Y, weights, offset, init = init,
                           controlvals = control, dist = dlist, scale = scale, nstrat = nstrata,
                           strata, parms = parms)
        break
      }
          
      init = init + runif(length(init), -1, 1)*init 
      n.init <- n.init + 1
     }
  }
  if (is.character(fit))
    fit <- list(fail = fit)
  else {
    if (scale == 0) {
      nvar <- length(fit$coefficients) - nstrata
      fit$scale <- exp(fit$coefficients[-(1:nvar)])
      if (nstrata == 1)
        names(fit$scale) <- NULL
      else names(fit$scale) <- levels(strata.keep)
      fit$coefficients <- fit$coefficients[1:nvar]
      fit$idf <- 1 + nstrata
    }
    else {
      fit$scale <- scale
      fit$idf <- 1
    }
    fit$loglik <- fit$loglik + logcorrect
  }
  if (!score)
    fit$score <- NULL
  fit$df.residual <- n - sum(fit$df)
  fit$terms <- Terms
  fit$contrasts <- contr.save
  if (length(xlevels))
    fit$xlevels <- xlevels
  fit$means <- apply(X, 2, mean)
  if (!is.null(weights))
    fit$weights <- weights
  fit$call <- Call
  fit$dist <- dist
  if (model)
    fit$model <- m
  if (x)
    fit$x <- X
  if (y)
    fit$y <- Ysave
  if (length(parms))
    fit$parms <- parms
  if (robust) {
    fit$naive.var <- fit$var
    if (!model)
      fit$model <- m
    if (length(cluster))
      fit$var <- crossprod(rowsum(residuals.survreg(fit,
                                                    "dfbeta"), cluster))
    else fit$var <- crossprod(residuals.survreg(fit, "dfbeta"))
    if (!model)
      fit$model <- NULL
  }
  singular <- (diag(fit$var) == 0)[1:length(fit$coefficients)]
  if (any(singular))
    fit$coefficients[singular] <- NA
  na.action <- attr(m, "na.action")
  if (length(na.action))
    fit$na.action <- na.action
  if (any(pterms))
    class(fit) <- c("survreg.penal", "survreg")
  else class(fit) <- "survreg"
  fit
}

survreg.fit <- function (x, y, weights, offset, init, controlvals, dist, scale = 0, 
          nstrat = 1, strata, parms = NULL, assign) 
{
  iter.max <- controlvals$iter.max
  eps <- controlvals$rel.tolerance
  toler.chol <- controlvals$toler.chol
  if (!is.matrix(x)) 
    stop("Invalid X matrix ")
  n <- nrow(x)
  nvar <- ncol(x)
  ny <- ncol(y)
  if (is.null(offset)) 
    offset <- rep(0, n)
  if (missing(weights) || is.null(weights)) 
    weights <- rep(1, n)
  else if (any(weights <= 0)) 
    stop("Invalid weights, must be >0")
  if (scale < 0) 
    stop("Invalid scale")
  if (scale > 0 && nstrat > 1) 
    stop("Cannot have both a fixed scale and strata")
  if (nstrat > 1 && (missing(strata) || length(strata) != n)) 
    stop("Invalid strata variable")
  if (nstrat == 1) 
    strata <- rep(1, n)
  if (scale > 0) 
    nstrat2 <- 0
  else nstrat2 <- nstrat
  if (is.character(dist)) {
    sd <- survreg.distributions[[dist]]
    if (is.null(sd)) 
      stop("Unrecognized distribution")
  }
  else sd <- dist
  if (!is.function(sd$density)) 
    stop("Missing density function in the definition of the distribution")
  dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
  if (is.na(dnum)) {
    dnum <- 4
    n2 <- n + sum(y[, ny] == 3)
    f.expr <- quote({
      if (length(parms)) temp <- sd$density(z, parms) else temp <- sd$density(z)
      if (!is.matrix(temp) || any(dim(temp) != c(n2, 5)) || 
          !is.numeric(temp)) stop("Density function returned an invalid matrix")
      as.vector(as.double(temp))
    })
    rho <- new.env()
  }
  else {
    f.expr <- 1
    rho <- 1
  }
  derfun <- function(y, eta, sigma, density, parms) {
    ny <- ncol(y)
    status <- y[, ny]
    z <- (y[, 1] - eta)/sigma
    dmat <- density(z, parms)
    dtemp <- dmat[, 3] * dmat[, 4]
    if (any(status == 3)) {
      z2 <- (y[, 2] - eta)/sigma
      dmat2 <- density(z2, parms)
    }
    else {
      dmat2 <- matrix(0, 1, 5)
      z2 <- 0
    }
    tdenom <- ((status == 0) * dmat[, 2]) + ((status == 1) * 
                                               1) + ((status == 2) * dmat[, 1]) + ((status == 3) * 
                                                                                     ifelse(z > 0, dmat[, 2] - dmat2[, 2], dmat2[, 1] - 
                                                                                              dmat[, 1]))
    tdenom <- 1/(tdenom * sigma)
    dg <- -tdenom * (((status == 0) * (0 - dmat[, 3])) + 
                       ((status == 1) * dmat[, 4]) + ((status == 2) * dmat[, 
                                                                           3]) + ((status == 3) * (dmat2[, 3] - dmat[, 3])))
    ddg <- (tdenom/sigma) * (((status == 0) * (0 - dtemp)) + 
                               ((status == 1) * dmat[, 5]) + ((status == 2) * dtemp) + 
                               ((status == 3) * (dmat2[, 3] * dmat2[, 4] - dtemp)))
    list(dg = dg, ddg = ddg - dg^2)
  }
  rescaled <- FALSE
  if (is.null(init) && all(x[, 1] == 1) && ncol(x) > 1) {
    okay <- apply(x, 2, function(z) all(z == 0 | z == 1))
    if (!all(okay)) {
      rescaled <- TRUE
      center <- ifelse(okay, 0, colMeans(x))
      stdev <- ifelse(okay, 1, apply(x, 2, sd))
      x <- scale(x, center, stdev)
    }
  }
  nvar2 <- nvar + nstrat2
  meanonly <- (nvar == 1 && all(x == 1))
  if (!meanonly) {
    yy <- ifelse(y[, ny] != 3, y[, 1], (y[, 1] + y[, 2])/2)
    coef <- sd$init(yy, weights, parms)
    if (scale > 0) 
      vars <- log(scale)
    else vars <- log(4 * coef[2])/2
    coef <- c(coef[1], rep(vars, nstrat))
    deriv <- derfun(y, yy, exp(vars), sd$density, parms)
    wt <- -1 * deriv$ddg * weights
    coef[1] <- sum(weights * deriv$dg + wt * (yy - offset))/sum(wt)
    fit0 <- .Call(survival:::Csurvreg6, iter = as.integer(20), nvar = as.integer(1), 
                  as.double(y), as.integer(ny), x = as.double(rep(1, 
                                                                  n)), as.double(weights), as.double(offset), coef = as.double(coef), 
                  as.integer(nstrat2), as.integer(strata), as.double(eps), 
                  as.double(toler.chol), as.integer(dnum), f.expr, 
                  rho)
  }
  if (is.numeric(init)) {
    if (length(init) == nvar && (nvar2 > nvar)) {
      init <- c(init, fit0$coef[-1])
    }
    if (length(init) != nvar2) 
      stop("Wrong length for initial parameters")
    if (scale > 0) 
      init <- c(init, log(scale))
  }
  else {
    if (meanonly) {
      yy <- ifelse(y[, ny] != 3, y[, 1], (y[, 1] + y[, 
                                                     2])/2)
      coef <- sd$init(yy, weights, parms)
      if (scale > 0) 
        vars <- rep(log(scale), nstrat)
      else vars <- rep(log(4 * coef[2])/2, nstrat)
    }
    else vars <- fit0$coef[-1]
    eta <- yy - offset
    deriv <- derfun(y, yy, exp(vars[strata]), sd$density, 
                    parms)
    wt <- -1 * deriv$ddg * weights
    coef <- coxph.wtest(t(x) %*% (wt * x), c((wt * eta + 
                                                weights * deriv$dg) %*% x), toler.chol = toler.chol)$solve
    init <- c(coef, vars)
  }
  fit <- .Call(survival:::Csurvreg6, iter = as.integer(iter.max), as.integer(nvar), 
               as.double(y), as.integer(ny), as.double(x), as.double(weights), 
               as.double(offset), as.double(init), as.integer(nstrat2), 
               as.integer(strata), as.double(eps), as.double(toler.chol), 
               as.integer(dnum), f.expr, rho)
  if (iter.max > 1 && fit$flag > nvar2) {
    # Commenting out so can run with warnings converted to errors without this
    # flagging
    # warning("Ran out of iterations and did not converge")
  }
  cname <- dimnames(x)[[2]]
  if (is.null(cname)) 
    cname <- paste("x", 1:ncol(x))
  if (scale == 0) 
    cname <- c(cname, rep("Log(scale)", nstrat))
  if (scale > 0) 
    fit$coef <- fit$coef[1:nvar2]
  names(fit$coef) <- cname
  if (meanonly) {
    coef0 <- fit$coef
    loglik <- rep(fit$loglik, 2)
  }
  else {
    coef0 <- fit0$coef
    names(coef0) <- c("Intercept", rep("Log(scale)", nstrat))
    loglik <- c(fit0$loglik, fit$loglik)
  }
  var.new <- matrix(fit$var, nvar2, dimnames = list(cname, 
                                                    cname))
  coef.new <- fit$coef
  if (rescaled) {
    vtemp <- diag(c(1/stdev, rep(1, nrow(var.new) - length(stdev))))
    vtemp[1, 2:nvar] <- -center[2:nvar]/stdev[2:nvar]
    coef.new <- drop(vtemp %*% coef.new)
    var.new <- vtemp %*% var.new %*% t(vtemp)
    names(coef.new) <- names(fit$coef)
  }
  temp <- list(coefficients = coef.new, icoef = coef0, var = var.new, 
               loglik = loglik, iter = fit$iter, linear.predictors = c(x %*% 
                                                                         fit$coef[1:nvar] + offset), df = length(fit$coef), 
               score = fit$u)
  temp
}

##############################


# Test cases ### Need to run through line 26 on plotting file

# # Test cases ### Need to run through line 26 on plotting file

# Scenario12 <- bind_rows(SubPlotDataList[33:48])
# Scenario22 <- bind_rows(SubPlotDataList[65:80])
# 
# # Single variable comparison
# ConfusedTable.Fun(Scenario12, Scenario22, "Distribution")
# ConfusedTable.Fun(Scenario12, Scenario22, "Regressors")
# 
# # Composite variable comparison
# ConfusedTable.Fun(Scenario12, Scenario22, "Regressors", "Distribution")
# 
# # Matching distributions only, then comparing regressors
# ConfusedTable.Fun(Scenario12, Scenario22, "Distribution", "Regressors", match_dist = TRUE)