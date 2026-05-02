### Simulating Survival Data
# Load packages and functions
source("./Survival Project/R Code/Survival - Functions and Packages.R")

# Initialize parameters
n <- 250
# c(Intercept, Treatment), c(Covariates)
Beta <- c(c(0.1, 0.2), c(0.15, 0.075, 0)) 

# Initialize settings
nsim <- 30000
Info.Criteria <- "AIC"
MS.Type <- "backward"

# Distribution list, censoring list, correlation list
SimulationResults <- ConfidenceIntervals <- Tabulated.Results <- Chosen.Mod.Results <- list()
index <- 0

# Run it
for(correlation in c("Independent", "ModerateCorrelation")){
  if (correlation == "Independent") {
    # Covariates are independent
    Sigma <- diag(length(Beta) - 2)
    
  } else if (correlation == "ModerateCorrelation") {
    Sigma <- matrix(2/3, nrow = length(Beta) - 2, ncol = length(Beta) - 2)
    diag(Sigma) <- 1
  }
  
  for(censoring in c("LowCens", "HighCens")){
    # Setting censoring level
    if (censoring == "LowCens") {
      # ~ 20 - 30% censored 
      CensPar <-  list(rate = 0.20)
      
    } else if (censoring == "HighCens") {
      # ~ 40 - 60% censored 
      CensPar <-  list(rate = 0.50)
    }
    
    for(distribution in c("weibull", "exponential", "lognormal", "cox")){
      # Setting distribution parameters
      if (distribution == "weibull") {
        ObservationParameters <- list(shape = 1.5, scale = 1.5)
        
      } else if (distribution == "exponential") {
        ObservationParameters <- list(rate = 0.50)
        
      } else if (distribution == "lognormal") {
        ObservationParameters <- list(shape = 0.25)
        
      } else if (distribution == "cox") {
        ObservationParameters <- list(shape = 0.50, scale = 0.50)
        
      } else {
        stop("Uh oh")
      }
      
      # Initialize storage objects
      IC.Object2 <- IC.Object <- matrix(NA, nrow = nsim, ncol = 4)
      Chosen.Model.Type2 <- Chosen.Dist2 <- Chosen.Model.Type <- Chosen.Dist <- list()
      Chosen.Model <- data.frame(matrix(NA, nrow = nsim, ncol = 2))
      
      # -------------------------------------
      # Full model objects
      # -------------------------------------
      Full.Cox.Coef <- Full.AFT.Coef <- Full.AFT2.Coef <- Full.AFT3.Coef <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Full.Cox.Coef.SE <- Full.AFT.Coef.SE <- Full.AFT2.Coef.SE <- Full.AFT3.Coef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Full.Cox.Z.Val <- Full.AFT.Z.Val <- Full.AFT2.Z.Val <- Full.AFT3.Z.Val <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Full.Cox.Concord <- matrix(NA, nrow = nsim, ncol = 2)
      Full.AFT.Scale <- Full.AFT3.Scale <- vector("double", length = nsim)
      Full.Cox.CI.Array <- Full.AFT.CI.Array <- Full.AFT2.CI.Array <- Full.AFT3.CI.Array <- array(NA, dim = c(length(Beta), 2, nsim))
      
      # -------------------------------------
      # Route 1 Part 1 objects (MS.XXX.XXX)
      # -------------------------------------
      MS.AFT.Coef <- MS.AFT2.Coef <- MS.AFT3.Coef <- MS.Cox.Coef <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.AFT.Coef.SE <- MS.AFT2.Coef.SE <- MS.AFT3.Coef.SE <- MS.Cox.Coef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.AFT.Z.Val <- MS.AFT2.Z.Val <- MS.AFT3.Z.Val <- MS.Cox.Z.Val <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.AFT.Scale <- MS.AFT3.Scale <- vector("double", length = nsim)
      MS.Cox.CI.Array <- MS.AFT.CI.Array <- MS.AFT2.CI.Array <- MS.AFT3.CI.Array <- array(NA, dim = c(length(Beta), 2, nsim))
      
      # -------------------------------------
      # Route 1 Part 2 objects (MS.Sel.XXX)
      # -------------------------------------
      MS.Sel.Cox.Coef <- MS.Sel.AFT.Coef <- MS.Sel.AFT2.Coef <- MS.Sel.AFT3.Coef <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.Sel.Cox.Coef.SE <- MS.Sel.AFT.Coef.SE <- MS.Sel.AFT2.Coef.SE <- MS.Sel.AFT3.Coef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.Sel.Cox.Z.Val <- MS.Sel.AFT.Z.Val <- MS.Sel.AFT2.Z.Val <- MS.Sel.AFT3.Z.Val <- matrix(NA, nrow = nsim, ncol = length(Beta))
      MS.Sel.AFT.Scale <- MS.Sel.AFT3.Scale <- vector("double", length = nsim)
      MS.Sel.Cox.CI.Array <- MS.Sel.AFT.CI.Array <- MS.Sel.AFT2.CI.Array <- MS.Sel.AFT3.CI.Array <- array(NA, dim = c(length(Beta), 2, nsim))
      
      # -------------------------------------
      # Route 2 Part 1 objects (Sel.XXX.XXX)
      # -------------------------------------
      Sel.Cox.Coef <- Sel.AFT.Coef <- Sel.AFT2.Coef <- Sel.AFT3.Coef <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.Cox.Coef.SE <- Sel.AFT.Coef.SE <- Sel.AFT2.Coef.SE <- Sel.AFT3.Coef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.Cox.Z.Val <- Sel.AFT.Z.Val <- Sel.AFT2.Z.Val <- Sel.AFT3.Z.Val <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.Cox.Concord <- matrix(NA, nrow = nsim, ncol = 2)
      Sel.AFT.Scale <- Sel.AFT3.Scale <- vector("double", length = nsim)
      Sel.Cox.CI.Array <- Sel.AFT.CI.Array <- Sel.AFT2.CI.Array <- Sel.AFT3.CI.Array <- array(NA, dim = c(length(Beta), 2, nsim))
      
      # -------------------------------------
      # Route 2 Part 2 objects (Sel.MS.XXX)
      # -------------------------------------
      Sel.MS.Cox.Coef <- Sel.MS.AFT.Coef <- Sel.MS.AFT2.Coef <- Sel.MS.AFT3.Coef <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.MS.Cox.Coef.SE <- Sel.MS.AFT.Coef.SE <- Sel.MS.AFT2.Coef.SE <- Sel.MS.AFT3.Coef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.MS.Cox.Z.Val <- Sel.MS.AFT.Z.Val <- Sel.MS.AFT2.Z.Val <- Sel.MS.AFT3.Z.Val <- matrix(NA, nrow = nsim, ncol = length(Beta))
      Sel.MS.AFT.Scale <- Sel.MS.AFT3.Scale <- vector("double", length = nsim)
      Sel.MS.Cox.CI.Array <- Sel.MS.AFT.CI.Array <- Sel.MS.AFT2.CI.Array <- Sel.MS.AFT3.CI.Array <- array(NA, dim = c(length(Beta), 2, nsim))
      
      
      for(i in 1:nsim) {
        tryCatch({
          # Set / change seed (starting from 1)
          set.seed(i)
          
          # Generate data
          GeneratedData <- Surv.Data.Gen.Fun(n = n, Beta = Beta, Sigma = Sigma, 
                                             ObsDist = distribution, ObsPar = ObservationParameters,
                                             CensDist = "exponential", CensPar = CensPar)
          
          # Fit Full AFT Models
          # AFT Model 1
          AFT.Mod <- survreg(Surv(time, delta, type = "right") ~ ., data = GeneratedData, 
                             dist = "weibull", beta = Beta,
                             control = survreg.control(maxiter = 500),
                             init = Beta)
          AFT.Full.Mod <- head(summary(AFT.Mod)[["table"]], -1)
          ind <- match(attr(AFT.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Full.AFT.Coef[i, c(1, ind + 1)] <- AFT.Full.Mod[, 1]
          Full.AFT.Coef.SE[i, c(1, ind + 1)] <- AFT.Full.Mod[, 2]
          Full.AFT.Z.Val[i, c(1, ind + 1)] <- AFT.Full.Mod[, 3]
          Full.AFT.Scale[i] <- AFT.Mod$scale
          Full.AFT.CI.Array[c(1, ind + 1), , i] <- confint(AFT.Mod) # CI Information
          
          # AFT Model 2
          AFT2.Mod <- survreg(Surv(time, delta, type = "right") ~ ., data = GeneratedData, 
                              dist = "exponential", beta = Beta,
                              control = survreg.control(maxiter = 500), 
                              init = Beta)
          AFT2.Full.Mod <- summary(AFT2.Mod)[["table"]]
          ind <- match(attr(AFT2.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Full.AFT2.Coef[i, c(1, ind + 1)] <- AFT2.Full.Mod[, 1]
          Full.AFT2.Coef.SE[i, c(1, ind + 1)] <- AFT2.Full.Mod[, 2]
          Full.AFT2.Z.Val[i, c(1, ind + 1)] <- AFT2.Full.Mod[, 3]
          Full.AFT2.CI.Array[c(1, ind + 1), , i] <- confint(AFT2.Mod) # CI Information
          
          # AFT Model 3
          AFT3.Mod <- survreg(Surv(time, delta, type = "right") ~ ., data = GeneratedData, 
                              dist = "lognormal", beta = Beta,
                              control = survreg.control(maxiter = 500), 
                              init = Beta)
          AFT3.Full.Mod <- head(summary(AFT3.Mod)[["table"]], -1)
          ind <- match(attr(AFT3.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Full.AFT3.Coef[i, c(1, ind + 1)] <- AFT3.Full.Mod[, 1]
          Full.AFT3.Coef.SE[i, c(1, ind + 1)] <- AFT3.Full.Mod[, 2]
          Full.AFT3.Z.Val[i, c(1, ind + 1)] <- AFT3.Full.Mod[, 3]
          Full.AFT3.Scale[i] <- AFT3.Mod$scale
          Full.AFT3.CI.Array[c(1, ind + 1), , i] <- confint(AFT3.Mod) # CI Information
          
          # Fit Full Cox Model
          Cox.Mod <- coxph(Surv(time, delta, type = "right") ~ ., data = GeneratedData)
          Cox.Full.Mod <- summary(Cox.Mod)
          ind <- match(attr(Cox.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1, 2)])
          
          # Saving relevant statistics
          Full.Cox.Coef[i, ind + 1] <- Cox.Full.Mod$coefficients[, 1]
          Full.Cox.Coef.SE[i, ind + 1] <- Cox.Full.Mod$coefficients[, 2]
          Full.Cox.Z.Val[i, ind + 1] <- Cox.Full.Mod$coefficients[, 3]
          Full.Cox.Concord[i, ] <- Cox.Full.Mod$concordance
          Full.Cox.CI.Array[ind + 1, , i] <- confint(Cox.Mod) # CI Information
          
          ### Route 1 - Pick which type of model then variable selection
          set.seed(i)
          
          ## Part 1 - Model Selection
          # Comparing the 4 models using AIC or BIC
          if(Info.Criteria == "AIC"){
            IC.Object[i, ] <- sapply(list(AFT.Mod = AFT.Mod, AFT2.Mod = AFT2.Mod,
                                          AFT3.Mod = AFT3.Mod, Cox.Mod = Cox.Mod), AIC)
          } else if (Info.Criteria == "BIC") {
            IC.Object[i, ] <- sapply(list(AFT.Mod = AFT.Mod, AFT2.Mod = AFT2.Mod,
                                          AFT3.Mod = AFT3.Mod, Cox.Mod = Cox.Mod), BIC)
          } else {
            stop("Information criteria ain't it bro")
          }
          
          # Which model survives to have variable selection done on it
          if (which.min(IC.Object[i, ]) == 1) {
            ChosenModel <- AFT.Mod
            Chosen.Dist[[i]] <- "weibull"
            Chosen.Model.Type[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(AFT.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            MS.AFT.Coef[i, c(1, ind + 1)] <- Full.AFT.Coef[i, c(1, ind + 1)]
            MS.AFT.Coef.SE[i, c(1, ind + 1)] <- Full.AFT.Coef.SE[i, c(1, ind + 1)]
            MS.AFT.Z.Val[i, c(1, ind + 1)] <- Full.AFT.Z.Val[i, c(1, ind + 1)]
            MS.AFT.Scale[i] <- Full.AFT.Scale[i]
            MS.AFT.CI.Array[c(1, ind + 1), , i] <- confint(ChosenModel) # CI Information
            
          } else if (which.min(IC.Object[i, ]) == 2) {
            ChosenModel <- AFT2.Mod
            Chosen.Dist[[i]] <- "exponential"
            Chosen.Model.Type[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(AFT2.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            MS.AFT2.Coef[i, c(1, ind + 1)] <- Full.AFT2.Coef[i, c(1, ind + 1)]
            MS.AFT2.Coef.SE[i, c(1, ind + 1)] <- Full.AFT2.Coef.SE[i, c(1, ind + 1)]
            MS.AFT2.Z.Val[i, c(1, ind + 1)] <- Full.AFT2.Z.Val[i, c(1, ind + 1)]
            MS.AFT2.CI.Array[c(1, ind + 1), , i] <- confint(ChosenModel) # CI Information
            
          } else if (which.min(IC.Object[i, ]) == 3) {
            ChosenModel <- AFT3.Mod
            Chosen.Dist[[i]] <- "lognormal"
            Chosen.Model.Type[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(AFT3.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            MS.AFT3.Coef[i, c(1, ind + 1)] <- Full.AFT3.Coef[i, c(1, ind + 1)]
            MS.AFT3.Coef.SE[i, c(1, ind + 1)] <- Full.AFT3.Coef.SE[i, c(1, ind + 1)]
            MS.AFT3.Z.Val[i, c(1, ind + 1)] <- Full.AFT3.Z.Val[i, c(1, ind + 1)]
            MS.AFT3.Scale[i] <- Full.AFT3.Scale[i]
            MS.AFT3.CI.Array[c(1, ind + 1), , i] <- confint(ChosenModel) # CI Information
            
          } else if (which.min(IC.Object[i, ]) == 4) {
            ChosenModel <- Cox.Mod
            Chosen.Dist[[i]] <- NA
            Chosen.Model.Type[[i]] <- "Cox"
            
            # Re-establish indices
            ind <- match(attr(Cox.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1, 2)])
            
            # Saving relevant statistics
            MS.Cox.Coef[i, ind + 1] <- Full.Cox.Coef[i, ind + 1]
            MS.Cox.Coef.SE[i, ind + 1] <- Full.Cox.Coef.SE[i, ind + 1]
            MS.Cox.Z.Val[i, ind + 1] <- Full.Cox.Z.Val[i, ind + 1]
            MS.Cox.CI.Array[ind + 1, , i] <- confint(ChosenModel) # CI Information
            
          } else {
            stop("Something has gone awry")
          }
          
          ## Variable Selection
          Sel.Mod <- Survival.MS.Fun(Data = GeneratedData, Beta = Beta, Model = Chosen.Model.Type[[i]], MS.Type = MS.Type, 
                                     IC = Info.Criteria, dist = Chosen.Dist[[i]])
          
          # Pulling stuff out from selected model
          if(attributes(Sel.Mod)[["class"]] == "coxph"){
            # Selected Model
            Cox.Sel.Mod <- summary(Sel.Mod)
            ind <- match(attr(Sel.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1, 2)])
            
            # Saving relevant statistics
            MS.Sel.Cox.Coef[i, ind + 1] <- Cox.Sel.Mod$coefficients[, 1]
            MS.Sel.Cox.Coef.SE[i, ind + 1] <- Cox.Sel.Mod$coefficients[, 2]
            MS.Sel.Cox.Z.Val[i, ind + 1] <- Cox.Sel.Mod$coefficients[, 3]
            MS.Sel.Cox.CI.Array[ind + 1, , i] <- confint(Sel.Mod) # CI Information
            
          } else if(attributes(Sel.Mod)[["class"]] == "survreg") {
            if(Chosen.Dist[[i]] == "weibull"){
              # Selected Model
              AFT.Sel.Mod <- head(summary(Sel.Mod)[["table"]], -1)
              ind <- match(attr(Sel.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
              
              # Saving relevant statistics
              MS.Sel.AFT.Coef[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 1]
              MS.Sel.AFT.Coef.SE[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 2]
              MS.Sel.AFT.Z.Val[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 3]
              MS.Sel.AFT.Scale[i] <- Sel.Mod$scale
              MS.Sel.AFT.CI.Array[c(1, ind + 1), , i] <- confint(Sel.Mod) # CI Information
              
            } else if (Chosen.Dist[[i]] == "exponential") {
              # Selected Model
              AFT.Sel.Mod <- summary(Sel.Mod)[["table"]]
              ind <- match(attr(Sel.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
              
              # Saving relevant statistics
              MS.Sel.AFT2.Coef[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 1]
              MS.Sel.AFT2.Coef.SE[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 2]
              MS.Sel.AFT2.Z.Val[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 3]
              MS.Sel.AFT2.CI.Array[c(1, ind + 1), , i] <- confint(Sel.Mod) # CI Information
              
            } else if(Chosen.Dist[[i]] == "lognormal") {
              # Selected Model
              AFT.Sel.Mod <- head(summary(Sel.Mod)[["table"]], -1)
              ind <- match(attr(Sel.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
              
              # Saving relevant statistics
              MS.Sel.AFT3.Coef[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 1]
              MS.Sel.AFT3.Coef.SE[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 2]
              MS.Sel.AFT3.Z.Val[i, c(1, ind + 1)] <- AFT.Sel.Mod[, 3]
              MS.Sel.AFT3.Scale[i] <- Sel.Mod$scale
              MS.Sel.AFT3.CI.Array[c(1, ind + 1), , i] <- confint(Sel.Mod) # CI Information
              
            } else {
              stop("cry")
            }
          } else {
            stop("The world is falling apart")
          }
          
          ### Route 2 - Variable selection then pick which type of model
          set.seed(i)
          # Perform AFT variable selections
          # AFT Model 1
          MS.AFT.Mod <- Survival.MS.Fun(Data = GeneratedData, Beta = Beta, Model = "AFT", MS.Type = MS.Type, 
                                        IC = Info.Criteria, dist = "weibull")
          MS.AFT.Mod.tbl <- head(summary(MS.AFT.Mod)[["table"]], -1)
          ind <- match(attr(MS.AFT.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Sel.AFT.Coef[i, c(1, ind + 1)] <- MS.AFT.Mod.tbl[, 1]
          Sel.AFT.Coef.SE[i, c(1, ind + 1)] <- MS.AFT.Mod.tbl[, 2]
          Sel.AFT.Z.Val[i, c(1, ind + 1)] <- MS.AFT.Mod.tbl[, 3]
          Sel.AFT.Scale[i] <- MS.AFT.Mod$scale
          Sel.AFT.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT.Mod) # CI Information
          
          # AFT Model 2
          MS.AFT2.Mod <- Survival.MS.Fun(Data = GeneratedData, Beta = Beta, Model = "AFT", MS.Type = MS.Type, 
                                         IC = Info.Criteria, dist = "exponential")
          MS.AFT2.Mod.tbl <- summary(MS.AFT2.Mod)[["table"]]
          ind <- match(attr(MS.AFT2.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Sel.AFT2.Coef[i, c(1, ind + 1)] <- MS.AFT2.Mod.tbl[, 1]
          Sel.AFT2.Coef.SE[i, c(1, ind + 1)] <- MS.AFT2.Mod.tbl[, 2]
          Sel.AFT2.Z.Val[i, c(1, ind + 1)] <- MS.AFT2.Mod.tbl[, 3]
          Sel.AFT2.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT2.Mod) # CI Information
          
          # AFT Model 3
          MS.AFT3.Mod <- Survival.MS.Fun(Data = GeneratedData, Beta = Beta, Model = "AFT", MS.Type = MS.Type, 
                                         IC = Info.Criteria, dist = "lognormal")
          MS.AFT3.Mod.tbl <- head(summary(MS.AFT3.Mod)[["table"]], -1)
          ind <- match(attr(MS.AFT3.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
          
          # Saving relevant statistics
          Sel.AFT3.Coef[i, c(1, ind + 1)] <- MS.AFT3.Mod.tbl[, 1]
          Sel.AFT3.Coef.SE[i, c(1, ind + 1)] <- MS.AFT3.Mod.tbl[, 2]
          Sel.AFT3.Z.Val[i, c(1, ind + 1)] <- MS.AFT3.Mod.tbl[, 3]
          Sel.AFT3.Scale[i] <- MS.AFT3.Mod$scale
          Sel.AFT3.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT3.Mod) # CI Information
          
          # Perform Cox variable selection
          MS.Cox.Mod <- Survival.MS.Fun(Data = GeneratedData, Beta = Beta, Model = "Cox", MS.Type = MS.Type, 
                                        IC = Info.Criteria)
          MS.Cox.Mod.Sum <- summary(MS.Cox.Mod)
          ind <- match(attr(MS.Cox.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1, 2)])
          
          # Saving relevant statistics
          Sel.Cox.Coef[i, ind + 1] <- MS.Cox.Mod.Sum$coefficients[, 1]
          Sel.Cox.Coef.SE[i, ind + 1] <- MS.Cox.Mod.Sum$coefficients[, 2]
          Sel.Cox.Z.Val[i, ind + 1] <- MS.Cox.Mod.Sum$coefficients[, 3]
          Sel.Cox.CI.Array[ind + 1, , i] <- confint(MS.Cox.Mod) # CI Information
          
          # Comparing the 4 models using AIC or BIC
          if(Info.Criteria == "AIC"){
            IC.Object2[i, ] <- sapply(list(MS.AFT.Mod = MS.AFT.Mod, MS.AFT2.Mod = MS.AFT2.Mod, 
                                           MS.AFT3.Mod = MS.AFT3.Mod, MS.Cox.Mod = MS.Cox.Mod), AIC)
          } else if (Info.Criteria == "BIC") {
            IC.Object2[i, ] <- sapply(list(MS.AFT.Mod = MS.AFT.Mod, MS.AFT2.Mod = MS.AFT2.Mod, 
                                           MS.AFT3.Mod = MS.AFT3.Mod, MS.Cox.Mod = MS.Cox.Mod), BIC)
          } else {
            stop("Information criteria not there man...")
          }
          
          # Which model is the best based on IC after model selection?
          if (which.min(IC.Object2[i, ]) == 1) {
            Chosen.Dist2[[i]] <- "weibull"
            Chosen.Model.Type2[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(MS.AFT.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            Sel.MS.AFT.Coef[i, c(1, ind + 1)] <- Sel.AFT.Coef[i, c(1, ind + 1)]
            Sel.MS.AFT.Coef.SE[i, c(1, ind + 1)] <- Sel.AFT.Coef.SE[i, c(1, ind + 1)]
            Sel.MS.AFT.Z.Val[i, c(1, ind + 1)] <- Sel.AFT.Z.Val[i, c(1, ind + 1)]
            Sel.MS.AFT.Scale[i] <- Sel.AFT.Scale[i]
            Sel.MS.AFT.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT.Mod) # CI Information
            
          } else if (which.min(IC.Object2[i, ]) == 2) {
            Chosen.Dist2[[i]] <- "exponential"
            Chosen.Model.Type2[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(MS.AFT2.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            Sel.MS.AFT2.Coef[i, c(1, ind + 1)] <- Sel.AFT2.Coef[i, c(1, ind + 1)] 
            Sel.MS.AFT2.Coef.SE[i, c(1, ind + 1)] <- Sel.AFT2.Coef.SE[i, c(1, ind + 1)] 
            Sel.MS.AFT2.Z.Val[i, c(1, ind + 1)] <- Sel.AFT2.Z.Val[i, c(1, ind + 1)]
            Sel.MS.AFT2.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT2.Mod) # CI Information
            
          } else if (which.min(IC.Object2[i, ]) == 3) {
            Chosen.Dist2[[i]] <- "lognormal"
            Chosen.Model.Type2[[i]] <- "AFT"
            
            # Re-establish indices
            ind <- match(attr(MS.AFT3.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1,2)])
            
            # Saving relevant statistics
            Sel.MS.AFT3.Coef[i, c(1, ind + 1)] <- Sel.AFT3.Coef[i, c(1, ind + 1)]
            Sel.MS.AFT3.Coef.SE[i, c(1, ind + 1)] <- Sel.AFT3.Coef.SE[i, c(1, ind + 1)]
            Sel.MS.AFT3.Z.Val[i, c(1, ind + 1)] <- Sel.AFT3.Z.Val[i, c(1, ind + 1)]
            Sel.MS.AFT3.Scale[i] <- Sel.AFT3.Scale[i]
            Sel.MS.AFT3.CI.Array[c(1, ind + 1), , i] <- confint(MS.AFT3.Mod) # CI Information
            
          } else if (which.min(IC.Object2[i, ]) == 4) {
            Chosen.Dist2[[i]] <- NA
            Chosen.Model.Type2[[i]] <- "Cox"
            
            # Re-establish indices
            ind <- match(attr(MS.Cox.Mod$terms, "term.labels"), colnames(GeneratedData)[-c(1, 2)])
            
            # Saving relevant statistics
            Sel.MS.Cox.Coef[i, ind + 1] <- Sel.Cox.Coef[i, ind + 1]
            Sel.MS.Cox.Coef.SE[i, ind + 1] <- Sel.Cox.Coef.SE[i, ind + 1]
            Sel.MS.Cox.Z.Val[i, ind + 1] <- Sel.Cox.Z.Val[i, ind + 1]
            Sel.MS.Cox.CI.Array[ind + 1, , i] <- confint(MS.Cox.Mod) # CI Information
            
          } else {
            stop("Something broke")
          }
          
          # Store for table of agreement
          Chosen.Model[i, ] <- c(Chosen.Dist[[i]], Chosen.Dist2[[i]]) 
          
          # Progress indicator
          if(i %% (nsim/5) == 0){
            print(i)
          }
          
        }, error = function(e) {
        message("Skipping iteration ", i, " due to error: ", e$message,
                " on simulation setting: ", correlation, "/", censoring, "/", distribution)
      })
      
      }
      
      # Tabulate chosen distributions and model types for both scenarios
      Chosen.Distributions <- table(unlist(Chosen.Dist))
      Chosen.Distributions2 <- table(unlist(Chosen.Dist2))
      Chosen.Model.Type <- table(unlist(Chosen.Model.Type))
      Chosen.Model.Type2 <- table(unlist(Chosen.Model.Type2))
      
      # Create a list organizing storage objects by model type
      Sims <- list(
        AFT1_Weibull = lapply(list(
          Full_Coef = Full.AFT.Coef,
          MS_Coef = MS.AFT.Coef,
          MS_Sel_Coef = MS.Sel.AFT.Coef,
          Sel_Coef = Sel.AFT.Coef,
          Sel_MS_Coef = Sel.MS.AFT.Coef,
          
          Full_Coef_SE = Full.AFT.Coef.SE,
          MS_Coef_SE = MS.AFT.Coef.SE,
          MS_Sel_Coef_SE = MS.Sel.AFT.Coef.SE,
          Sel_Coef_SE = Sel.AFT.Coef.SE,
          Sel_MS_Coef_SE = Sel.MS.AFT.Coef.SE,
          
          Full_Z_Val = Full.AFT.Z.Val,
          MS_Z_Val = MS.AFT.Z.Val,
          MS_Sel_Z_Val = MS.Sel.AFT.Z.Val,
          Sel_Z_Val = Sel.AFT.Z.Val,
          Sel_MS_Z_Val = Sel.MS.AFT.Z.Val,
          
          Full_Scale = Full.AFT.Scale,
          MS_Scale = MS.AFT.Scale,
          MS_Sel_Scale = MS.Sel.AFT.Scale,
          Sel_Scale = Sel.AFT.Scale,
          Sel_MS_Scale = Sel.MS.AFT.Scale
        ), ApplyColnames),
        
        AFT2_Exponential = lapply(list(
          Full_Coef = Full.AFT2.Coef,
          MS_Coef = MS.AFT2.Coef,
          MS_Sel_Coef = MS.Sel.AFT2.Coef,
          Sel_Coef = Sel.AFT2.Coef,
          Sel_MS_Coef = Sel.MS.AFT2.Coef,
          
          Full_Coef_SE = Full.AFT2.Coef.SE,
          MS_Coef_SE = MS.AFT2.Coef.SE,
          MS_Sel_Coef_SE = MS.Sel.AFT2.Coef.SE,
          Sel_Coef_SE = Sel.AFT2.Coef.SE,
          Sel_MS_Coef_SE = Sel.MS.AFT2.Coef.SE,
          
          Full_Z_Val = Full.AFT2.Z.Val,
          MS_Z_Val = MS.AFT2.Z.Val,
          MS_Sel_Z_Val = MS.Sel.AFT2.Z.Val,
          Sel_Z_Val = Sel.AFT2.Z.Val,
          Sel_MS_Z_Val = Sel.MS.AFT2.Z.Val
        ), ApplyColnames),
        
        AFT3_Lognormal = lapply(list(
          Full_Coef = Full.AFT3.Coef,
          MS_Coef = MS.AFT3.Coef,
          MS_Sel_Coef = MS.Sel.AFT3.Coef,
          Sel_Coef = Sel.AFT3.Coef,
          Sel_MS_Coef = Sel.MS.AFT3.Coef,
          
          Full_Coef_SE = Full.AFT3.Coef.SE,
          MS_Coef_SE = MS.AFT3.Coef.SE,
          MS_Sel_Coef_SE = MS.Sel.AFT3.Coef.SE,
          Sel_Coef_SE = Sel.AFT3.Coef.SE,
          Sel_MS_Coef_SE = Sel.MS.AFT3.Coef.SE,
          
          Full_Z_Val = Full.AFT3.Z.Val,
          MS_Z_Val = MS.AFT3.Z.Val,
          MS_Sel_Z_Val = MS.Sel.AFT3.Z.Val,
          Sel_Z_Val = Sel.AFT3.Z.Val,
          Sel_MS_Z_Val = Sel.MS.AFT3.Z.Val,
          
          Full_Scale = Full.AFT3.Scale,
          MS_Scale = MS.AFT3.Scale,
          MS_Sel_Scale = MS.Sel.AFT3.Scale,
          Sel_Scale = Sel.AFT3.Scale,
          Sel_MS_Scale = Sel.MS.AFT3.Scale
        ), ApplyColnames),
        
        Cox_Regression = lapply(list(
          Full_Coef = Full.Cox.Coef,
          MS_Coef = MS.Cox.Coef,
          MS_Sel_Coef = MS.Sel.Cox.Coef,
          Sel_Coef = Sel.Cox.Coef,
          Sel_MS_Coef = Sel.MS.Cox.Coef,
          
          Full_Coef_SE = Full.Cox.Coef.SE,
          MS_Coef_SE = MS.Cox.Coef.SE,
          MS_Sel_Coef_SE = MS.Sel.Cox.Coef.SE,
          Sel_Coef_SE = Sel.Cox.Coef.SE,
          Sel_MS_Coef_SE = Sel.MS.Cox.Coef.SE,
          
          Full_Z_Val = Full.Cox.Z.Val,
          MS_Z_Val = MS.Cox.Z.Val,
          MS_Sel_Z_Val = MS.Sel.Cox.Z.Val,
          Sel_Z_Val = Sel.Cox.Z.Val,
          Sel_MS_Z_Val = Sel.MS.Cox.Z.Val,
          
          Full_Concord = Full.Cox.Concord,
          Sel_Concord = Sel.Cox.Concord
        ), ApplyColnames),
        
        Miscellaneous = lapply(list(
          IC_Object = IC.Object,
          IC_Object2 = IC.Object2,
          Chosen_Model_Type = Chosen.Model.Type,
          Chosen_Model_Type2 = Chosen.Model.Type2,
          Chosen_Dist = Chosen.Distributions,
          Chosen_Dist2 = Chosen.Distributions2,
          Chosen.Model = Chosen.Model
        ), ApplyColnames)
      )
      
      # Confidence Interval Arrays
      CI.Sims <- list(
        AFT1_Weibull = list(
          Full_CI = Full.AFT.CI.Array,
          MS_CI = MS.AFT.CI.Array,
          MS_Sel_CI = MS.Sel.AFT.CI.Array,
          Sel_CI = Sel.AFT.CI.Array,
          Sel_MS_CI = Sel.MS.AFT.CI.Array
        ),
        
        AFT2_Exponential = list(
          Full_CI = Full.AFT2.CI.Array,
          MS_CI = MS.AFT2.CI.Array,
          MS_Sel_CI = MS.Sel.AFT2.CI.Array,
          Sel_CI = Sel.AFT2.CI.Array,
          Sel_MS_CI = Sel.MS.AFT2.CI.Array
        ),
        
        AFT3_Lognormal = list(
          Full_CI = Full.AFT3.CI.Array,
          MS_CI = MS.AFT3.CI.Array,
          MS_Sel_CI = MS.Sel.AFT3.CI.Array,
          Sel_CI = Sel.AFT3.CI.Array,
          Sel_MS_CI = Sel.MS.AFT3.CI.Array
        ),
        
        Cox_Regression = list(
          Full_CI = Full.Cox.CI.Array,
          MS_CI = MS.Cox.CI.Array,
          MS_Sel_CI = MS.Sel.Cox.CI.Array,
          Sel_CI = Sel.Cox.CI.Array,
          Sel_MS_CI = Sel.MS.Cox.CI.Array
        )
      )
      
      # Performing coverage calculations on the fly
      ConfidenceIntervals[[correlation]][[censoring]][[distribution]] <- CI.Sims
      
      # Store all the crap for current simulation settings
      SimulationResults[[correlation]][[censoring]][[distribution]] <- Sims
      
      # Aggregate results for current simulation settings
      index <- index + 1
      Tabulated.Results[[index]] <- SimulationTabulation.Fun(SimResults = SimulationResults,
                                                             CI = ConfidenceIntervals,
                                                             correlation = correlation,
                                                             censoring = censoring, 
                                                             distribution = distribution)[[1]]
      Chosen.Mod.Results[[index]] <- SimulationTabulation.Fun(SimResults = SimulationResults,
                                                             CI = ConfidenceIntervals,
                                                             correlation = correlation,
                                                             censoring = censoring, 
                                                             distribution = distribution)[[2]]
    }
  }
}

# Monster output
# SimulationResults
# Tabulated.Results
 
# Save simulation results to file
# save(SimulationResults, file = "./Survival Project/Survival Simulation Results/SurvSimResultsRData")
# save(Tabulated.Results, file = "./Survival Project/Survival Simulation Results/SurvSimTabulated.RData")
# save(Chosen.Mod.Results, file = "./Survival Project/Survival Simulation Results/SurvSimChosenMod.RData")

