### Mixture Distribution Simulations
set.seed(10262025)
library(mvtnorm)
library(tidyverse)
library(glmnet)
library(gridExtra)
library(RColorBrewer)

# Function generating sampling distributions post-model selection
Mixed.Fun <- function(mu.vec = c(2, 4, 6),                                  # Mean vector for covariates
                      cov.mat = rbind(c(5, 4, 5), c(4, 6, 5), c(5, 5, 7)),  # Covariance matrix
                      n = 200,                                              # Sample Size
                      epsilon.error = 10,                                   # Error variance
                      beta.seq = c(3, 0, 1, 2),                             # Coefficients (including intercept)
                      MS.type = "forward",                                  # Model selection type (forward, backwards, step-wise)
                      nsim = 10000,                                         # Number of simulations
                      criteria = "AIC",                                     # Model selection criterion (AIC or BIC)
                      plot.type = "Density",                                # Automatic plot type (Density, Histogram, 2D*)
                      plot.stat = "T3",                                     # Which statistic to automatically plot
                      correct.model.plot = FALSE,                           # Whether to include correct plot or not
                      L.R. = FALSE,                                         # To perform LASSO and Ridge or not
                      data.split = FALSE,                                   # To perform Data splitting or not
                      SP = 0.70,                                            # Data splitting ratio (training data proportion)
                      verbose = FALSE){                                     # Print progress if desired
  # Error Check
  if(data.split && L.R.){
   stop("Data-Splitting is not currently supported for LASSO and Ridge") 
  }
  
  ### Empty object creation
  Reg.Sel.List <- rep(NA, nsim)
  
  # Objects for model selected model
  res.SE <- rep(NA, nsim)
  Coef.mat <- CoefSE.mat <- Tval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  
  # Objects for "correct" model
  Cres.SE <- rep(NA, nsim)
  CCoef.mat <- CCoefSE.mat <- CTval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  
  # Objects for full model
  Fres.SE <- rep(NA, nsim)
  FCoef.mat <- FCoefSE.mat <- FTval.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
  
  # Objects for LASSO and Ridge Models
  LCoef.mat <- RCoef.mat <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
  
  if(data.split){
    ## Initialize objects for train and test data
    Reg.Sel.List.train <- rep(NA, nsim)
    Reg.Sel.List.test <- rep(NA, nsim)
    
    # Objects for model selected model for train and test
    res.SE.train <- rep(NA, nsim)
    res.SE.test <- rep(NA, nsim)
    
    Coef.mat.train <- CoefSE.mat.train <- Tval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    Coef.mat.test <- CoefSE.mat.test <- Tval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    
    # Objects for "correct" model for train and test
    Cres.SE.train <- rep(NA, nsim)
    Cres.SE.test <- rep(NA, nsim)
    
    CCoef.mat.train <- CCoefSE.mat.train <- CTval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    CCoef.mat.test <- CCoefSE.mat.test <- CTval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    
    # Objects for full model for train and test
    Fres.SE.train <- rep(NA, nsim)
    Fres.SE.test <- rep(NA, nsim)
    
    FCoef.mat.train <- FCoefSE.mat.train <- FTval.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    FCoef.mat.test <- FCoefSE.mat.test <- FTval.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq))
    
    # Objects for LASSO and Ridge Models for train and test
    LCoef.mat.train <- RCoef.mat.train <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
    LCoef.mat.test <- RCoef.mat.test <- matrix(data = NA, nrow = nsim, ncol = length(beta.seq) - 1)
  }
  
  ### Start Simulations
  for(i in 1:nsim){
  ## Random MVN Sample
  x.mat <- rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
  
  # ABC labeling of regressors unless there is a lot
  if(length(mu.vec) > 24){
    colnames(x.mat) <- paste("X", 1:length(mu.vec), sep = "")
  }
  colnames(x.mat) <- LETTERS[1:length(mu.vec)]
  
  # Creating response
  y.vec <- beta.seq[1] + beta.seq[-1]%*%t(x.mat) + rnorm(n, mean = 0, sd = sqrt(epsilon.error))

  # Creating dataframe
  da.data <- cbind.data.frame(Y = as.vector(y.vec), x.mat)
  
  # Performing specified model selection technique
  if(MS.type %in% c("forward", "backward", "both")){
  null.model <- lm(Y ~ 1, data = da.data)
  Mod.Select <- step(object = null.model, scope = as.formula(paste("~", paste(colnames(x.mat), collapse = " + "))), 
                     direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
  
  # Noting which regressors were chosen in the selected mode
  Reg.Sel.List[i] <- paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = "")
  
  # Pulling summary from selected model
  lm.select <- summary(lm(formula = Mod.Select$call$formula, data = da.data))
  
  # Checking which index the variables chosen correspond to
  ind <- match(attr(Mod.Select$terms, "term.labels"), colnames(x.mat))
  
  # Saving coefficients, coefficient standard errors, t-statistics and model's residual SE
  Coef.mat[i, c(1, ind+1)] <- lm.select$coefficients[,1]
  CoefSE.mat[i, c(1, ind+1)] <- lm.select$coefficients[,2]
  Tval.mat[i, c(1, ind+1)] <- lm.select$coefficients[,3]
  res.SE[i] <- lm.select$sigma
  
  # Pulling summary from"correct" model
  # Checking for which regressors should be in the model
  Correct.Regs <- colnames(x.mat)[which(beta.seq[-1] != 0)]
  lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = da.data))
  indd <- match(attr(lm.correct$terms, "term.labels"), colnames(x.mat))

  # Saving relevant statistics
  CCoef.mat[i, c(1, indd+1)] <- lm.correct$coefficients[,1]
  CCoefSE.mat[i, c(1, indd+1)] <- lm.correct$coefficients[,2]
  CTval.mat[i, c(1, indd+1)] <- lm.correct$coefficients[,3]
  Cres.SE[i] <- lm.correct$sigma

  # Pulling summary from full model
  lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = da.data))
  
  # Saving relevant statistics
  FCoef.mat[i, 1:length(beta.seq)] <- lm.full$coefficients[,1]
  FCoefSE.mat[i, 1:length(beta.seq)] <- lm.full$coefficients[,2]
  FTval.mat[i, 1:length(beta.seq)] <- lm.full$coefficients[,3]
  Fres.SE[i] <- lm.full$sigma
  }
  
  # Data-Splitting Strategy - Performing specified model selection technique
  if(data.split == TRUE){
    # Split data into training and test sets
    index <- sample(1:n, size = n * SP)
    train.data <- da.data[index, ]
    test.data <- da.data[-index, ]
    
    # Perform model selection only on training data
    if(MS.type %in% c("forward", "backward", "both")) {
      ## Train Data
      # Null model for model selection (on training data)
      null.model <- lm(Y ~ 1, data = train.data)
      Mod.Select <- step(object = null.model, scope = as.formula(paste("~", paste(colnames(x.mat), collapse = " + "))), 
                         direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
      
      # Noting which regressors were chosen in the selected model (for train data)
      Reg.Sel.List.train[i] <- paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = "")
      
      # Pulling summary from selected model (on train data)
      lm.select <- summary(lm(formula = Mod.Select$call$formula, data = train.data))
      
      # Checking which index the variables chosen correspond to
      ind <- match(attr(Mod.Select$terms, "term.labels"), colnames(x.mat))
      
      # Saving coefficients, coefficient standard errors, t-statistics, and model's residual SE (for train data)
      Coef.mat.train[i, c(1, ind+1)] <- lm.select$coefficients[,1]
      CoefSE.mat.train[i, c(1, ind+1)] <- lm.select$coefficients[,2]
      Tval.mat.train[i, c(1, ind+1)] <- lm.select$coefficients[,3]
      res.SE.train[i] <- lm.select$sigma
      
      # Pulling summary from "correct" model (on train data)
      Correct.Regs <- colnames(x.mat)[which(beta.seq[-1] != 0)]
      lm.correct <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = train.data))
      indd <- match(attr(lm.correct$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the correct model (for train data)
      CCoef.mat.train[i, c(1, indd+1)] <- lm.correct$coefficients[,1]
      CCoefSE.mat.train[i, c(1, indd+1)] <- lm.correct$coefficients[,2]
      CTval.mat.train[i, c(1, indd+1)] <- lm.correct$coefficients[,3]
      Cres.SE.train[i] <- lm.correct$sigma
      
      # Pulling summary from full model (on train data)
      lm.full <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = train.data))
      
      # Saving relevant statistics for the full model (for train data)
      FCoef.mat.train[i, 1:length(beta.seq)] <- lm.full$coefficients[,1]
      FCoefSE.mat.train[i, 1:length(beta.seq)] <- lm.full$coefficients[,2]
      FTval.mat.train[i, 1:length(beta.seq)] <- lm.full$coefficients[,3]
      Fres.SE.train[i] <- lm.full$sigma
      
      ## Test Data
      # Apply the selected model to the test data
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
      Coef.mat.test[i, c(1, ind.test+1)] <- lm.select.test$coefficients[,1]
      CoefSE.mat.test[i, c(1, ind.test+1)] <- lm.select.test$coefficients[,2]
      Tval.mat.test[i, c(1, ind.test+1)] <- lm.select.test$coefficients[,3]
      res.SE.test[i] <- lm.select.test$sigma
      
      # Pulling summary from "correct" model (on test data)
      lm.correct.test <- summary(lm(formula = as.formula(paste("Y ~", paste(Correct.Regs, collapse = " + "))), data = test.data))
      indd.test <- match(attr(lm.correct.test$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics for the correct model (for test data)
      CCoef.mat.test[i, c(1, indd.test+1)] <- lm.correct.test$coefficients[,1]
      CCoefSE.mat.test[i, c(1, indd.test+1)] <- lm.correct.test$coefficients[,2]
      CTval.mat.test[i, c(1, indd.test+1)] <- lm.correct.test$coefficients[,3]
      Cres.SE.test[i] <- lm.correct.test$sigma
      
      # Pulling summary from full model (on test data)
      lm.full.test <- summary(lm(formula = as.formula(paste("Y ~", paste(colnames(x.mat), collapse = " + "))), data = test.data))
      
      # Saving relevant statistics for the full model (for test data)
      FCoef.mat.test[i, 1:length(beta.seq)] <- lm.full.test$coefficients[,1]
      FCoefSE.mat.test[i, 1:length(beta.seq)] <- lm.full.test$coefficients[,2]
      FTval.mat.test[i, 1:length(beta.seq)] <- lm.full.test$coefficients[,3]
      Fres.SE.test[i] <- lm.full.test$sigma
    }
  }
  
  if(L.R. == TRUE){
  # LASSO 
  # Cross validation to calculate lambda
  LASSO.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 1, nlambda = 100, standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  LASSO.mod <- glmnet(x = x.mat, y = y.vec, alpha = 1, lambda = LASSO.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  LCoef.mat[i, 1:(length(beta.seq)-1)] <- as.matrix(LASSO.mod$beta)
  
  # Ridge
  # Cross validation to calculate lambda
  Ridge.CV <- cv.glmnet(x = x.mat, y = y.vec, alpha = 0, nlambda = 100, standardize = TRUE, thresh = 1e-10)
  # Fitting Model
  Ridge.mod <- glmnet(x = x.mat, y = y.vec, alpha = 0, lambda = Ridge.CV$lambda.min, 
                      standardize = TRUE, thresh = 1e-10)
  # Pulling coefficients
  RCoef.mat[i, 1:(length(beta.seq)-1)] <- as.matrix(Ridge.mod$beta)
  }
  
  # Print progress if desired
  if (verbose) {
    if(i %% (nsim / 10) == 0){
    cat(i, "Simulations completed", "\n")
    }
  }
  }
      
  ### Aggregating Results
  # Combining results for post-model selection values
  Results <- cbind.data.frame(Reg.Sel.List, Coef.mat, CoefSE.mat, Tval.mat, res.SE)
  
  # Results for "correct" model
  CResults <- cbind.data.frame(Reg.Sel.List, CCoef.mat, CCoefSE.mat, CTval.mat, Cres.SE)
  
  # Results for full model
  FResults <- cbind.data.frame(Reg.Sel.List, FCoef.mat, FCoefSE.mat, FTval.mat, Fres.SE)
  
  # Results for LASSO and Ridge Models
  LRResults <- cbind.data.frame(LCoef.mat, RCoef.mat)
  
  # Labeling
  colnames(LRResults) <- c(paste("LBeta", 1:(length(beta.seq)-1), sep = ""),
                           paste("RBeta", 1:(length(beta.seq)-1), sep = ""))
  colnames(Results) <- c("Regressors", paste("Beta", 0:(length(beta.seq)-1), sep = ""), 
                         paste("SE.Beta", 0:(length(beta.seq)-1), sep = ""),
                         paste("T", 0:(length(beta.seq)-1), sep = ""), "Res.SE")
  colnames(FResults) <- colnames(CResults) <- colnames(Results) # [-1]
  
  if(data.split){
    # Aggregating Results for Train and Test
    # Combine results for the model selected from the training data
    Train.Results <- cbind.data.frame(Reg.Sel.List.train, Coef.mat.train, CoefSE.mat.train, Tval.mat.train, res.SE.train)
    Test.Results <- cbind.data.frame(Reg.Sel.List.test, Coef.mat.test, CoefSE.mat.test, Tval.mat.test, res.SE.test)
    
    # Combine results for the "correct" model
    Train.CResults <- cbind.data.frame(CCoef.mat.train, CCoefSE.mat.train, CTval.mat.train, Cres.SE.train)
    Test.CResults <- cbind.data.frame(CCoef.mat.test, CCoefSE.mat.test, CTval.mat.test, Cres.SE.test)
    
    # Combine results for the full model
    Train.FResults <- cbind.data.frame(FCoef.mat.train, FCoefSE.mat.train, FTval.mat.train, Fres.SE.train)
    Test.FResults <- cbind.data.frame(FCoef.mat.test, FCoefSE.mat.test, FTval.mat.test, Fres.SE.test)
    
    # Labeling for Train and Test Results
    colnames(Train.Results) <- c("Regressors", paste("Beta", 0:(length(beta.seq)-1), sep = ""), 
                                 paste("SE.Beta", 0:(length(beta.seq)-1), sep = ""),
                                 paste("T", 0:(length(beta.seq)-1), sep = ""), "Res.SE")
    colnames(Test.Results) <- colnames(Train.Results)
    
    # Labeling for Correct Model Results (Train and Test)
    colnames(Train.CResults) <- colnames(Train.Results)[-1]
    colnames(Test.CResults) <- colnames(Train.Results)[-1]
    
    # Labeling for Full Model Results (Train and Test)
    colnames(Train.FResults) <- colnames(Train.Results)[-1]
    colnames(Test.FResults) <- colnames(Train.Results)[-1]
  }
  
  ### Automatic Plotting
  if (plot.type == "Density") {
    # Calculate the range of all layers
    x_values <- c(
      Results[[plot.stat]],
      FResults[[plot.stat]],
      if (correct.model.plot && !all(is.na(CResults[[plot.stat]]))) CResults[[plot.stat]] else NULL,
      if (startsWith(plot.stat, "Beta")) beta.seq[as.numeric(sub("Beta", "", plot.stat)) + 1] else NULL
    )
    x_range <- range(x_values, na.rm = TRUE)
    
    # Expand the range slightly
    x_range <- c(x_range[1] - 0.075 * diff(x_range), x_range[2] + 0.075 * diff(x_range))
    
    # Prepare LaTeX formatted title and axis labels
    plot_title <- if (startsWith(plot.stat, "Beta")) {
      bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
    } else {
      paste("Density plot of", plot.stat)
    }
    
    plot_x_label <- if (startsWith(plot.stat, "Beta")) {
      bquote(hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
    } else {
      plot.stat
    }
    
    # Start building the plot
    plot <- ggplot() +
      geom_density(aes(x = Results[[plot.stat]], fill = "Selected Model"), 
                   color = "#5ab9e8", alpha = 0.10) +
      geom_density(aes(x = FResults[[plot.stat]], fill = "Full Model"), 
                   color = "#7beb78", alpha = 0.10) +
      labs(title = plot_title,
           subtitle = "Comparison of Selected Model with Full Model",
           x = plot_x_label, y = "Density") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center") +
      scale_fill_manual(values = c("Selected Model" = "#5ab9e8", "Full Model" = "#7beb78")) +
      guides(fill = guide_legend(title = "Sampling Distributions")) +
      scale_x_continuous(limits = x_range)
    
    # Add true coefficient vertical line
    if (startsWith(plot.stat, "Beta")) {
      beta_number <- as.numeric(sub("Beta", "", plot.stat))
      if (!is.na(beta_number)) {
        plot <- plot + 
          geom_vline(xintercept = beta.seq[beta_number + 1], linetype = "dashed")
      } else {
        warning("Invalid plot.stat format: Unable to extract Beta coefficient number.")
      }
    }
    
    # Add correct model layer if enabled
    if (correct.model.plot) {
      if (all(is.na(CResults[[plot.stat]]))) {
        message("Cannot plot correct model sampling distribution due to inexistence")
      } else {
        plot <- plot +
          geom_density(aes(x = CResults[[plot.stat]], fill = "Correct Model"), 
                       color = "#d55ae8", alpha = 0.10) +
          labs(subtitle = "Comparison of Selected, Full, and Correct Model") +
          scale_fill_manual(values = c(
            "Selected Model" = "#5ab9e8", 
            "Full Model" = "#7beb78", 
            "Correct Model" = "#d55ae8"
          ))
      }
    }
    
    # Print the plot
    print(plot)
    
  } else if(plot.type == "Histogram"){
    # Start with the base plot for histograms
    plot <- ggplot() +
      geom_histogram(aes(x = Results[[plot.stat]], fill = "Selected Model"), color = "#5ab9e8", 
                     bins = round(log(nsim)*10), alpha = 0.10) +
      geom_histogram(aes(x = FResults[[plot.stat]], fill = "Full Model"), color = "#7beb78", 
                     bins = round(log(nsim)*10), alpha = 0.10) +
      labs(title = paste("Histogram of", plot.stat),
           subtitle = "Comparison of Selected Model with Full Model",
           x = plot.stat,
           y = "Count") +
      theme_minimal() + 
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center") +
      scale_fill_manual(values = c("Selected Model" = "#5ab9e8", "Full Model" = "#7beb78")) +
      guides(fill = guide_legend(title = "Sampling Distributions"))
    
    # Conditionally add the layer for CResults (histogram)
    if (correct.model.plot == TRUE){
      plot <- plot + 
        geom_histogram(aes(x = CResults[[plot.stat]], fill = "Correct Model"), color = "#d55ae8", 
                       bins = round(log(nsim)*10), alpha = 0.10) + 
        labs(subtitle = "Comparison of Selected, Full, and Correct Model") +
        scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                     "Full Model" = "#7beb78", 
                                     "Correct Model" = "#d55ae8"))
    }
    # Print the plot
    print(plot)
  } else {
    # 2D Plot - Takes a vector of length 2 for the overhead function input 
    # - Probably not going to use this much right now, little effort here
    ggplot() +
      geom_density_2d_filled(aes(x = Results[[plot.stat[1]]], y = Results[[plot.stat[2]]]), alpha = 0.5) +
      labs(title = paste("Bivariate Plot of", plot.stat[1], "vs.", plot.stat[2]),
           x = plot.stat[1],
           y = plot.stat[2]) +
      theme_minimal() +
      theme(legend.position = "none")
    
  }
  # Output
  if(data.split){
    return(Output = list(
      Post.Model.Selection = Results, 
      Correct.Model = CResults, 
      Full.Model = FResults,
      Train.Post.Model.Selection = Train.Results, 
      Test.Post.Model.Selection = Test.Results,
      Train.Correct.Model = Train.CResults, 
      Test.Correct.Model = Test.CResults, 
      Train.Full.Model = Train.FResults, 
      Test.Full.Model = Test.FResults))
  } else{
    if(L.R. == TRUE){
      return(Output = list(Post.Model.Selection = Results, Correct.Model = CResults, 
                           Full.Model = FResults, LASSO.Ridge = LRResults))
    } else{
      return(Output = list(Post.Model.Selection = Results, Correct.Model = CResults, 
                           Full.Model = FResults))
    }
  }
}

### Some external plotting functions because there's a lot of plotting ahead...
# Plots of coefficients and their t-stats
Plot.Beta.Fun <- function(data, plot.stat = "Beta1") {
  # Prepare LaTeX formatted title and axis labels
  plot_title <- if (startsWith(plot.stat, "Beta")) {
    bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))] ~ "Sampling Distribution")
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (startsWith(plot.stat, "Beta")) {
    bquote(hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
  } else {
    plot.stat
  }
  
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[["Post.Model.Selection"]][plot.stat]), fill = "Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["Full.Model"]][plot.stat]), fill = "Full Model"), 
                 color = "#7beb78", alpha = 0.10)
  
  # Check if the data for the "Correct Model" is all NAs
  if (any(!is.na(data[["Correct.Model"]][plot.stat]))) {
    plot <- plot + geom_density(aes(x = as.matrix(data[["Correct.Model"]][plot.stat]), fill = "Correct Model"), 
                                color = "#d55ae8", alpha = 0.10)
  }
  
  # Add labels and theme to the plot
  plot <- plot + 
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                 "Full Model" = "#7beb78", 
                                 "Correct Model" = "#d55ae8")) +
    guides(fill = guide_legend(title = ""))
  
  # Add vertical line for Beta coefficient
  if (startsWith(plot.stat, "Beta")) {
    beta_number <- as.numeric(sub("Beta", "", plot.stat))
    if (!is.na(beta_number)) {
      plot <- plot + 
        geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed")
    } else {
      warning("Invalid plot.stat format: Unable to extract Beta coefficient number.")
    }
  }
  
  # Print the plot
  print(plot)
}
Plot.t.stat.Fun <- function(data, plot.stat = "T1") {
  # Prepare LaTeX formatted title and axis labels
  T.Num <- as.numeric(sub("T", "", plot.stat))
  plot_title <- if (startsWith(plot.stat, "T")) {
    bquote("Sampling Distribution Density Plot: t-Statistic of " ~ hat(beta)[.(T.Num)])
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (startsWith(plot.stat, "T")) {
    bquote("t-Statistic of " ~ hat(beta)[.(T.Num)])
  } else {
    plot.stat
  }
  
  # Initialize the plot
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[["Post.Model.Selection"]][plot.stat]), fill = "Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["Full.Model"]][plot.stat]), fill = "Full Model"), 
                 color = "#7beb78", alpha = 0.10)
  
  # Check if the data for the "Correct Model" is all NAs
  if (any(!is.na(data[["Correct.Model"]][plot.stat]))) {
    plot <- plot + geom_density(aes(x = as.matrix(data[["Correct.Model"]][plot.stat]), fill = "Correct Model"), 
                                color = "#d55ae8", alpha = 0.10)
  }
  
  # Add labels and theme to the plot
  plot <- plot + 
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                 "Full Model" = "#7beb78", 
                                 "Correct Model" = "#d55ae8")) +
    guides(fill = guide_legend(title = ""))
  
  # Print the plot
  print(plot)
}
# LASSO and Ridge Coefficients
Plot.LR.Beta.Fun <- function(data, plot.stat = "Beta1") {
  # Prepare LaTeX formatted title and axis labels
  plot_title <- if (grepl("Beta", plot.stat)) {
    bquote("Density plot of " ~ hat(beta)[.(as.numeric(gsub("\\D", "", plot.stat)))] ~ "Sampling Distribution")
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (grepl("Beta", plot.stat)) {
    bquote(hat(beta)[.(as.numeric(gsub("\\D", "", plot.stat)))])
  } else {
    plot.stat
  }
  
  # Create the plot
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[["Post.Model.Selection"]][plot.stat]), fill = "Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["Correct.Model"]][plot.stat]), fill = "Correct Model"), 
                 color = "#d55ae8", alpha = 0.10) + 
    geom_density(aes(x = as.matrix(data[["LASSO.Ridge"]][paste0("L", plot.stat)]), fill = "LASSO Model"), 
                 color = "#3dff7b", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["LASSO.Ridge"]][paste0("R", plot.stat)]), fill = "Ridge Model"), 
                 color = "#ff983d", alpha = 0.10) +
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                 "Correct Model" = "#d55ae8", 
                                 "LASSO Model" = "#3dff7b",
                                 "Ridge Model" = "#ff983d")) +
    guides(fill = guide_legend(title = ""))
  
  # Add vertical dashed line for Beta if applicable
  if (startsWith(plot.stat, "Beta")) {
    beta_number <- as.numeric(gsub("\\D", "", plot.stat)) # Extract beta number
    if (!is.na(beta_number) && beta_number <= length(Beta.seq)) { # Ensure Beta.seq is long enough
      plot <- plot + 
        geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed") # Directly use Beta.seq from the global environment
    } else {
      warning("Invalid beta_number or Beta.seq is too short.")
    }
  }
  
  print(plot)
}
# Plotting of train or test sampling distributions (not of use)
Plot.TT.Beta.Fun <- function(data, data.type = "Test", plot.stat = "Beta1") {
  # Prepare LaTeX formatted title and axis labels
  plot_title <- if (startsWith(plot.stat, "Beta")) {
    bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))] ~ "Sampling Distribution")
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (startsWith(plot.stat, "Beta")) {
    bquote(hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
  } else {
    plot.stat
  }
  
  # Dynamically create the names for data subsets based on the data.type (Train/Test)
  post_model_selection_col <- paste0(data.type, ".", "Post.Model.Selection")
  full_model_col <- paste0(data.type, ".", "Full.Model")
  correct_model_col <- paste0(data.type, ".", "Correct.Model")
  
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[[post_model_selection_col]][plot.stat]), fill = "Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[[full_model_col]][plot.stat]), fill = "Full Model"), 
                 color = "#7beb78", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[[correct_model_col]][plot.stat]), fill = "Correct Model"), 
                 color = "#d55ae8", alpha = 0.10) + 
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                 "Full Model" = "#7beb78", 
                                 "Correct Model" = "#d55ae8")) +
    guides(fill = guide_legend(title = ""))
  
  if (startsWith(plot.stat, "Beta")) {
    beta_number <- as.numeric(sub("Beta", "", plot.stat))
    if (!is.na(beta_number)) {
      plot <- plot + 
        geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed")
    } else {
      warning("Invalid plot.stat format: Unable to extract Beta coefficient number.")
    }
  }
  
  print(plot)
}
# Plotting train vs test (showing how much it's begin corrected)
Plot.DSTT.Beta.Fun <- function(data, plot.stat = "Beta1") {
  # Prepare LaTeX formatted title and axis labels
  plot_title <- if (startsWith(plot.stat, "Beta")) {
    bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))] ~ "Sampling Distribution")
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (startsWith(plot.stat, "Beta")) {
    bquote(hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
  } else {
    plot.stat
  }
  
  # Create the plot
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[["Train.Post.Model.Selection"]][plot.stat]), fill = "Data-Split Train Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["Correct.Model"]][plot.stat]), fill = "Correct Model"), 
                 color = "#d55ae8", alpha = 0.10) + 
    geom_density(aes(x = as.matrix(data[["Test.Post.Model.Selection"]][plot.stat]), fill = "Data-Split Test Model"), 
                 color = "#ff983d", alpha = 0.10) +
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Data-Split Train Selected Model" = "#5ab9e8", 
                                 "Correct Model" = "#d55ae8", 
                                 "Data-Split Test Model" = "#ff983d")) +
    guides(fill = guide_legend(title = ""))
  
  # Add vertical dashed line for Beta if applicable
  if (startsWith(plot.stat, "Beta")) {
    beta_number <- as.numeric(sub("Beta", "", plot.stat))
    if (!is.na(beta_number)) {
      plot <- plot + 
        geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed")
    } else {
      warning("Invalid plot.stat format: Unable to extract Beta coefficient number.")
    }
  }
  
  print(plot)
}
# Plotting the comparison of selected model with (test data) and without (full sim number of with mod. sel.)
# data splitting (and correct model for reference) (what we're interested in)
Plot.DS.Beta.Fun <- function(data, plot.stat = "Beta1") {
  # Prepare LaTeX formatted title and axis labels
  plot_title <- if (startsWith(plot.stat, "Beta")) {
    bquote("Density plot of " ~ hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))] ~ "Sampling Distribution")
  } else {
    paste("Density plot of", plot.stat)
  }
  
  plot_x_label <- if (startsWith(plot.stat, "Beta")) {
    bquote(hat(beta)[.(as.numeric(sub("Beta", "", plot.stat)))])
  } else {
    plot.stat
  }
  
  # Create the plot
  plot <- ggplot() +
    geom_density(aes(x = as.matrix(data[["Post.Model.Selection"]][plot.stat]), fill = "Selected Model"), 
                 color = "#5ab9e8", alpha = 0.10) +
    geom_density(aes(x = as.matrix(data[["Correct.Model"]][plot.stat]), fill = "Correct Model"), 
                 color = "#d55ae8", alpha = 0.10) + 
    geom_density(aes(x = as.matrix(data[["Test.Post.Model.Selection"]][plot.stat]), fill = "Data-Split Model"), 
                 color = "#ff983d", alpha = 0.10) +
    labs(title = plot_title,
         y = "Density",
         x = plot_x_label) +
    theme_minimal() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 1, color = "gray50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center") +
    scale_fill_manual(values = c("Selected Model" = "#5ab9e8", 
                                 "Correct Model" = "#d55ae8", 
                                 "Data-Split Model" = "#ff983d")) +
    guides(fill = guide_legend(title = ""))
  
  # Add vertical dashed line for Beta if applicable
  if (startsWith(plot.stat, "Beta")) {
    beta_number <- as.numeric(sub("Beta", "", plot.stat))
    if (!is.na(beta_number)) {
      plot <- plot + 
        geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed")
    } else {
      warning("Invalid plot.stat format: Unable to extract Beta coefficient number.")
    }
  }
  
  print(plot)
}
# Plotting the coefficient sampling distribution for each model from post-model selection
Plot.Decomp.Beta.Fun <- function(data, plot.stat = "Beta1", ncol = 2) {
  # Extract unique covariates from the Regressors column
  Mods <- names(table(data$Regressors))
  
  # Subset the data for each covariate
  subsets <- lapply(Mods, function(mod) {
    subset(data, Regressors == mod)
  })
  
  # Generate a color palette for the `fill`
  num_colors <- length(Mods)
  fill_colors <- brewer.pal(min(num_colors, 12), "Set3")
  
  # Create individual plots for each subset
  individual_plots <- lapply(seq_along(subsets), function(i) {
    if (any(!is.na(subsets[[i]][[plot.stat]]))) {
      plot <- ggplot(subsets[[i]], aes(x = subsets[[i]][[plot.stat]], fill = Mods[i])) +
        geom_density(alpha = 0.33, color = "black") +
        scale_fill_manual(values = fill_colors[i], name = "Regressor") +
        labs(title = paste("Model", Mods[i], "Density Plot", plot.stat, "Estimates"),
             x = "",
             y = "Density") +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "none"
        )
      
      # Add vertical line for Beta coefficient if Beta.seq is provided
      if (!is.null(Beta.seq) && startsWith(plot.stat, "Beta")) {
        beta_number <- as.numeric(sub("Beta", "", plot.stat))
        if (!is.na(beta_number) && beta_number + 1 <= length(Beta.seq)) {
          plot <- plot +
            geom_vline(xintercept = Beta.seq[beta_number + 1], linetype = "dashed", color = "black")
        } else {
          warning("Invalid Beta coefficient index or Beta.seq is too short.")
        }
      }
      
      return(plot)
    } else {
      NULL
    }
  })
  
  # Filter out any NULL plots
  individual_plots <- Filter(Negate(is.null), individual_plots)
  
  # Combine all individual plots
  combined_plot <- grid.arrange(grobs = individual_plots, ncol = ncol)
  
  # Return the combined plot
  return(combined_plot)
}
# Plotting the coefficient sampling distribution for each model from real model
Plot.Decomp.Beta.Fun.Correct <- function(data, Beta.seq, ncol = 2, colors = NULL) {
  # Validate inputs
  if (is.null(Beta.seq) || length(Beta.seq) < 4) {
    stop("Beta.seq must be provided and have at least 4 elements.")
  }
  
  # Set default colors if none are provided
  if (is.null(colors)) {
    colors <- c("#FF5733", "#33FF57", "#3357FF", "#F333FF") 
  }
  
  # Ensure enough colors are provided
  if (length(colors) < length(Beta.seq)) {
    stop("Not enough colors provided for the number of Beta coefficients.")
  }
  
  # Getting plotting dataframe
  subsets <- list(data[["Correct.Model"]]["Beta0"],
                  data[["Correc.Model"]]["Beta1"],
                  data[["Correc.Model"]]["Beta2"],
                  data[["Correc.Model"]]["Beta3"])
  
  # Create individual plots for each subset
  individual_plots <- lapply(seq_along(subsets), function(i) {
    current_data <- subsets[[i]]
    
    # Check if the subset has valid data
    if (any(!is.na(current_data[[1]]))) {  # Use the first column in each subset
      plot <- ggplot(current_data, aes(x = current_data[[1]])) +  # Reference the first column dynamically
        geom_density(alpha = 0.33, color = "black", fill = colors[i]) +
        labs(
          title = bquote("Density Plot of" ~ beta[.(i - 1)] ~ "Estimates"),
          x = "",
          y = "Density"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "none"
        )
      
      # Add vertical line for the current Beta coefficient
      if (i <= length(Beta.seq)) {
        plot <- plot +
          geom_vline(xintercept = Beta.seq[i], linetype = "dashed", color = "black")
      } else {
        warning("Beta.seq does not have enough elements for the current plot.")
      }
      
      return(plot)
    } else {
      NULL
    }
  })
  
  # Filter out any NULL plots
  individual_plots <- Filter(Negate(is.null), individual_plots)
  
  # Combine all individual plots
  combined_plot <- grid.arrange(grobs = individual_plots, ncol = ncol)
  
  # Return the combined plot
  return(combined_plot)
}

### Putting Reg.Sel.List at the beginning of the correct and full model results and the naming of columns
### Noting in case it fucks up something major in the future...

### Exploration 1
mean.vec <- c(1, 1, 1)
Sigma.Mat <- matrix(0.50, nrow = 3, ncol = 3)
diag(Sigma.Mat) <- 1
n <- 50
epsilon.var <- 16
Beta.seq <- c(2, -0.5, 0, 0.5)

# P1 - Generic Setting
Test.Fun1 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE)

# Breakdown of how often each model was chosen
table(Test.Fun1$Post.Model.Selection$Regressors)
Plot.Decomp.Beta.Fun(data = Test.Fun1$Post.Model.Selection, plot.stat = "Beta1", ncol = 2)
Plot.Decomp.Beta.Fun(data = Test.Fun1$Post.Model.Selection, plot.stat = "Beta2", ncol = 2)
Plot.Decomp.Beta.Fun(data = Test.Fun1$Post.Model.Selection, plot.stat = "Beta3", ncol = 2)

# Trying to breakdown mixture without model selection 
# Plot.Decomp.Beta.Fun(data = Test.Fun1$Correct.Model, plot.stat = "Beta1", ncol = 2)
# This is wrong, can't have beta1 estimates for model without beta 1, should fix later

Plot.Decomp.Beta.Fun.Correct(data = Test.Fun1, Beta.seq = Beta.seq, ncol = 2)

# Plotting mixtures of above
p1 <- Plot.Beta.Fun(data = Test.Fun1, plot.stat = "Beta1")
p2 <- Plot.Beta.Fun(data = Test.Fun1, plot.stat = "Beta2")
p3 <- Plot.Beta.Fun(data = Test.Fun1, plot.stat = "Beta3")
p4 <- Plot.t.stat.Fun(data = Test.Fun1, plot.stat = "T1")
p5 <- Plot.t.stat.Fun(data = Test.Fun1, plot.stat = "T2")
p6 <- Plot.t.stat.Fun(data = Test.Fun1, plot.stat = "T3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# P2 - Way up the sample size
n <- 300
Test.Fun1.2 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE)

p1 <- Plot.Beta.Fun(data = Test.Fun1.2, plot.stat = "Beta1")
p2 <- Plot.Beta.Fun(data = Test.Fun1.2, plot.stat = "Beta2")
p3 <- Plot.Beta.Fun(data = Test.Fun1.2, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, nrow = 1)

p4 <- Plot.t.stat.Fun(data = Test.Fun1.2, plot.stat = "T1")
p5 <- Plot.t.stat.Fun(data = Test.Fun1.2, plot.stat = "T2")
p6 <- Plot.t.stat.Fun(data = Test.Fun1.2, plot.stat = "T3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


# P3 - Xi uncorrelated
n <- 50
Sigma.Mat <- diag(1, nrow = 3, ncol = 3)
Test.Fun1.3 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                         epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                         MS.type = "forward", nsim = 10000, criteria = "AIC",
                         correct.model.plot = TRUE)

p1 <- Plot.Beta.Fun(data = Test.Fun1.3, plot.stat = "Beta1")
p2 <- Plot.Beta.Fun(data = Test.Fun1.3, plot.stat = "Beta2")
p3 <- Plot.Beta.Fun(data = Test.Fun1.3, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, nrow = 1)

p4 <- Plot.t.stat.Fun(data = Test.Fun1.3, plot.stat = "T1")
p5 <- Plot.t.stat.Fun(data = Test.Fun1.3, plot.stat = "T2")
p6 <- Plot.t.stat.Fun(data = Test.Fun1.3, plot.stat = "T3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


# P4 - Changing error variance
Sigma.Mat <- matrix(0.50, nrow = 3, ncol = 3)
diag(Sigma.Mat) <- 1
epsilon.var <- 36

Test.Fun1.4 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                         epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                         MS.type = "forward", nsim = 10000, criteria = "AIC",
                         correct.model.plot = TRUE)

p1 <- Plot.Beta.Fun(data = Test.Fun1.4, plot.stat = "Beta1")
p2 <- Plot.Beta.Fun(data = Test.Fun1.4, plot.stat = "Beta2")
p3 <- Plot.Beta.Fun(data = Test.Fun1.4, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, nrow = 1)

p4 <- Plot.t.stat.Fun(data = Test.Fun1.4, plot.stat = "T1")
p5 <- Plot.t.stat.Fun(data = Test.Fun1.4, plot.stat = "T2")
p6 <- Plot.t.stat.Fun(data = Test.Fun1.4, plot.stat = "T3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


### Exploration 2
mean.vec <- rep(1, 5)
Sigma.Mat <- matrix(0.50, nrow = 5, ncol = 5)
diag(Sigma.Mat) <- 1
n <- 50
epsilon.var <- 16
Beta.seq <- c(2, 0.0625, 0.125, 0.25, 0.50, 1.0)

Test.Fun2 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE)

p0 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T0")
p1 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T1")
p2 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T2")
p3 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T3")
p4 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T4")
p5 <- Plot.t.stat.Fun(data = Test.Fun2, plot.stat = "T5")
grid.arrange(p0, p1, p2, p3, p4, p5, nrow = 2)

# Again but center coefficients around 0 
Beta.seq <- c(2, -0.25, -0.125, 0, -0.125, 0.25)

Test.Fun2.2 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE)

p0 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T0")
p1 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T1")
p2 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T2")
p3 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T3")
p4 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T4")
p5 <- Plot.t.stat.Fun(data = Test.Fun2.2, plot.stat = "T5")
grid.arrange(p0, p1, p2, p3, p4, p5, nrow = 2)

### Exploration 3
mean.vec <- rep(1, 6)
Sigma.Mat <- matrix(0.50, nrow = 6, ncol = 6)
diag(Sigma.Mat) <- 1
n <- 50
epsilon.var <- 25
Beta.seq <- c(2, 0.125, 0.25, 0.50, 1.0, 2.0, 3.0)
Test.Fun3 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE, L.R. = TRUE)
p1 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta1")
p2 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta2")
p3 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta3")
p4 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta4")
p5 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta5")
p6 <- Plot.LR.Beta.Fun(data = Test.Fun3, plot.stat = "Beta6")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

### Exploration 4
# Was going to see if I could do this for a random effect sampling distribution... But my 
# GLMM function is not working ATM

### Exploration 5
set.seed(10262025)
mean.vec <- c(1, 1, 1)
Sigma.Mat <- matrix(0.5, nrow = 3, ncol = 3)
diag(Sigma.Mat) <- 1
n <- 50
epsilon.var <- 25
Beta.seq <- c(1, -3, -1, -2)

Test.Fun5 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE, SP = 0.5, data.split = TRUE)

# Something odd is occurring...
p1 <- Plot.DS.Beta.Fun(data = Test.Fun5, plot.stat = "Beta1")
p2 <- Plot.DS.Beta.Fun(data = Test.Fun5, plot.stat = "Beta2")
p3 <- Plot.DS.Beta.Fun(data = Test.Fun5, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, nrow = 1)

p4 <- Plot.DSTT.Beta.Fun(data = Test.Fun5, plot.stat = "Beta1")
p5 <- Plot.DSTT.Beta.Fun(data = Test.Fun5, plot.stat = "Beta2")
p6 <- Plot.DSTT.Beta.Fun(data = Test.Fun5, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# Changing proportion for data splitting
set.seed(10262025)
Test.Fun5.2 <- Mixed.Fun(mu.vec = mean.vec, cov.mat = Sigma.Mat, n = n, 
                       epsilon.error = epsilon.var, beta.seq = Beta.seq, 
                       MS.type = "forward", nsim = 10000, criteria = "AIC",
                       correct.model.plot = TRUE, SP = 0.75, data.split = TRUE)

p1 <- Plot.DS.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta1")
p2 <- Plot.DS.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta2")
p3 <- Plot.DS.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, nrow = 1)

p4 <- Plot.DSTT.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta1")
p5 <- Plot.DSTT.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta2")
p6 <- Plot.DSTT.Beta.Fun(data = Test.Fun5.2, plot.stat = "Beta3")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


### Function can't actually do backwards selection bro