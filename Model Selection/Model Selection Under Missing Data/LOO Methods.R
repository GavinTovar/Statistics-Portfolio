### Leave one out method evaluations
# Source functions and packages
source("./Missing Data Project/R Code/(ST 599) Functions and Packages.R")

### Setting
nsim <- 5
n <- 200
B <- 100
Beta <- c(0, c(0, 0.5, 1, 1.5, 2)) # 0 Intercept, no X1

### Initialize empty objects
# Scenario 1
NOIM.Boot.ChosenRegressors <- Stack.NOIM.Boot.MS <- NOIM.Jack.ChosenRegressors <- Stack.NO.MI.Jack.MS <- list()
# Initialize empty matrices for stability selection elements
StabSelectForwardCoef <- StabSelectForwardCoef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
StabSelectBestSubsetCoef <- StabSelectBestSubsetCoef.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))

# Scenario 2
Indiv.NOBS.MS <- list()
# Scenario 3
Std.MS <- list()
# Scenario 4
Chosen.Boot.Regressors <- Stack.Boot.MS <- Chosen.Jack.Regressors <- Stack.Jack.MS <- list()
# Initialize empty matrices for stability selection elements
ForwardSelection.KNN.Det.Imp.Data.Coef <- ForwardSelection.KNN.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
ForwardSelection.KNN.Stoc.Imp.Data.Coef <- ForwardSelection.KNN.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
ForwardSelection.Reg.Det.Imp.Data.Coef <- ForwardSelection.Reg.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
ForwardSelection.Reg.Stoc.Imp.Data.Coef <- ForwardSelection.Reg.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
ForwardSelection.Mean.Det.Imp.Data.Coef <- ForwardSelection.Mean.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
ForwardSelection.Mean.Stoc.Imp.Data.Coef <- ForwardSelection.Mean.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.KNN.Det.Imp.Data.Coef <- BestSubset.KNN.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.KNN.Stoc.Imp.Data.Coef <- BestSubset.KNN.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.Reg.Det.Imp.Data.Coef <- BestSubset.Reg.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.Reg.Stoc.Imp.Data.Coef <- BestSubset.Reg.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.Mean.Det.Imp.Data.Coef <- BestSubset.Mean.Det.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))
BestSubset.Mean.Stoc.Imp.Data.Coef <- BestSubset.Mean.Stoc.Imp.Data.SE <- matrix(NA, nrow = nsim, ncol = length(Beta))

### These objects should be initialized within the simulation setting loops below instead ################

# List for each missing mechanism simulation run
SimulationResults <- MissingPropResults <- Sims <- list()

for(correlation in c("Independent", "ModerateCorrelation")){
  if (correlation == "Independent") {
    # Covariates are independent
    Sigma <- diag(length(Beta) - 1)
    
  } else if (correlation == "ModerateCorrelation") {
    Sigma <- matrix(2/3, nrow = length(Beta) - 1, ncol = length(Beta) - 1)
    diag(Sigma) <- 1
  }
  
  for(MissingProportion in c("Low", "High")){
    # Setting missing proportion
    if (MissingProportion == "Low") {
      MissProportion <- 0.10
      
    } else if (MissingProportion == "High") {
      MissProportion <- 0.25
    }
    
    for(MissMechanism in c("MCAR", "MNAR", "MAR")){
      for(l in 1:nsim){
        set.seed(l)
        ### Generate data (Using the same for all situations so to be comparable)
        GeneratedData <- GenData.Fun(n = n, Beta = Beta, Sigma = Sigma, 
                                     MissMech = MissMechanism, MissProp = MissProportion)
        
        ##### 
        ### Situation 1 (No imputation)
        set.seed(l)
        ### Resampling Procedure
        # Bootstrap
        BootstrappedData <- Bootstrap.Fun(Data = GeneratedData, B = B)
        
        # Complete Case Cleaning
        CCBSData <- ArrayCleaning(BootstrappedData)
        
        ### Model Selection Section
        ## Bootstrap Section
        # Individual Model Selection
        Indiv.Boot.CC.MS <- MS.CC.Fun(Data = CCBSData, Beta = Beta, criteria = "AIC")
        
        # Stability Selection 
        NOIM.Boot.ChosenRegressors[[l]] <- Stability.Selection.CC.Fun(Indiv.Boot.CC.MS, pip = 0.50)
        
        # Fit model on generated data with selected variables
        Forward <- FitStabilitySelectedModel(Response = "Y", 
                                             Regressors = NOIM.Boot.ChosenRegressors[[l]]$Selection.Regs$ForwardSelection, 
                                             data = GeneratedData)
        # Saving coefficient values and their standard errors
        StabSelectForwardCoef[l, ] <- Forward[["coefficients"]]
        StabSelectForwardCoef.SE[l, ] <- Forward[["standard_errors"]]
        
        # BestSubsets Selection
        BestSubset <- FitStabilitySelectedModel(Response = "Y", 
                                             Regressors = NOIM.Boot.ChosenRegressors[[l]]$Selection.Regs$BestSubset, 
                                             data = GeneratedData)
        
        # Saving coefficient values and their standard errors
        StabSelectBestSubsetCoef[l, ] <- BestSubset[["coefficients"]]
        StabSelectBestSubsetCoef.SE[l, ] <- BestSubset[["standard_errors"]]
        
        # Stacking Model Selection
        Stack.NOIM.Boot.MS[[l]] <- MS.CC.Stack.Fun(Data = CCBSData, Beta = Beta, criteria = "AIC")
        
        #####
        ### Situation 2 (No resampling)
        set.seed(l)
        ### Imputation
        ImputedData <- Imp.NOBS.Fun(Data = GeneratedData, KNNsize = 3)
        
        ### Model Selection Section
        Indiv.NOBS.MS[[l]] <- MS.NOBS.Fun(Data = ImputedData, Beta = Beta, criteria = "AIC")
        
        # Notice: No stability selection or stacking needed
        
        #####
        ### Situation 3 (No resampling or imputation)
        set.seed(l)
        # Complete Case Cleaning
        Clean.Data <- GeneratedData[complete.cases(GeneratedData), ]
        
        ### Model Selection Section
        Std.MS[[l]] <- Std.MS.Fun(Data = Clean.Data, Beta = Beta, criteria = "AIC")
        
        #####
        ### Situation 4 (Resampling and imputation)
        set.seed(l)
        ### Resampling Procedure
        # Bootstrap
        BootstrappedData <- Bootstrap.Fun(Data = GeneratedData, B = B)
        
        ### Imputation
        # Bootstrap Imputation
        ImputedBootstraps <- Imp.Parallel.Fun(Data = BootstrappedData, KNNsize = 3)
        
        ### Model Selection Section
        ## Bootstrap Section
        # Individual Model Selection
        Indiv.Boot.MS <- MS.Fun(Data = ImputedBootstraps, Beta = Beta, criteria = "AIC")
        
        # Stability Selection 
        Chosen.Boot.Regressors[[l]] <- Stability.Selection.Fun(Indiv.Boot.MS, pip = 0.50)
        
        # Creating models of stability selected models
        StabSelections <- ApplySSModel(Selection.Regs = Chosen.Boot.Regressors[[l]]$Selection.Regs,
                                       Response = "Y",
                                       data = GeneratedData)
        
        # Storing coefficient results
        ForwardSelection.KNN.Det.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$KNN.Det.Imp.Data$coefficients
        ForwardSelection.KNN.Det.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$KNN.Det.Imp.Data$standard_errors
        
        ForwardSelection.KNN.Stoc.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$KNN.Stoc.Imp.Data$coefficients
        ForwardSelection.KNN.Stoc.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$KNN.Stoc.Imp.Data$standard_errors
        
        ForwardSelection.Reg.Det.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$Reg.Det.Imp.Data$coefficients
        ForwardSelection.Reg.Det.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$Reg.Det.Imp.Data$standard_errors
        
        ForwardSelection.Reg.Stoc.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$Reg.Stoc.Imp.Data$coefficients
        ForwardSelection.Reg.Stoc.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$Reg.Stoc.Imp.Data$standard_errors
        
        ForwardSelection.Mean.Det.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$Mean.Det.Imp.Data$coefficients
        ForwardSelection.Mean.Det.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$Mean.Det.Imp.Data$standard_errors
        
        ForwardSelection.Mean.Stoc.Imp.Data.Coef[l, ] <- StabSelections$ForwardSelection$Mean.Stoc.Imp.Data$coefficients
        ForwardSelection.Mean.Stoc.Imp.Data.SE[l, ] <- StabSelections$ForwardSelection$Mean.Stoc.Imp.Data$standard_errors
        
        BestSubset.KNN.Det.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$KNN.Det.Imp.Data$coefficients
        BestSubset.KNN.Det.Imp.Data.SE[l, ] <- StabSelections$BestSubset$KNN.Det.Imp.Data$standard_errors
        
        BestSubset.KNN.Stoc.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$KNN.Stoc.Imp.Data$coefficients
        BestSubset.KNN.Stoc.Imp.Data.SE[l, ] <- StabSelections$BestSubset$KNN.Stoc.Imp.Data$standard_errors
        
        BestSubset.Reg.Det.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$Reg.Det.Imp.Data$coefficients
        BestSubset.Reg.Det.Imp.Data.SE[l, ] <- StabSelections$BestSubset$Reg.Det.Imp.Data$standard_errors
        
        BestSubset.Reg.Stoc.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$Reg.Stoc.Imp.Data$coefficients
        BestSubset.Reg.Stoc.Imp.Data.SE[l, ] <- StabSelections$BestSubset$Reg.Stoc.Imp.Data$standard_errors
        
        BestSubset.Mean.Det.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$Mean.Det.Imp.Data$coefficients
        BestSubset.Mean.Det.Imp.Data.SE[l, ] <- StabSelections$BestSubset$Mean.Det.Imp.Data$standard_errors
        
        BestSubset.Mean.Stoc.Imp.Data.Coef[l, ] <- StabSelections$BestSubset$Mean.Stoc.Imp.Data$coefficients
        BestSubset.Mean.Stoc.Imp.Data.SE[l, ] <- StabSelections$BestSubset$Mean.Stoc.Imp.Data$standard_errors
        
        
        # Stacking Model Selection
        Stack.Boot.MS[[l]] <- MS.Stack.Fun(Data = ImputedBootstraps, Beta = Beta, criteria = "AIC")
        
        # Progress Statement
        if(l %% 5 == 0){
          cat(MissMechanism, "Simulations, iteration", l, "\n")
        }
      }
      
      # Collecting outputs
      SimResults <- list(
        Scenario1 = list(No.Imp.SS.Boot.Regs = NOIM.Boot.ChosenRegressors,
                         No.Imp.Stack.Boot.Regs = Stack.NOIM.Boot.MS,
                         StabilitySelectionForwardCoefs = StabSelectForwardCoef,
                         StabilitySelectionForwardCoefSEs = StabSelectForwardCoef.SE,
                         StabilitySelectionBestSubsetCoefs = StabSelectBestSubsetCoef,
                         StabilitySelectionBestSubsetCoefSEs = StabSelectBestSubsetCoef.SE),
        Scenario2 = list(Indiv.NOBS.MS),
        Scenario3 = list(Std.MS),
        Scenario4 = list(Chosen.Boot.Regressors = Chosen.Boot.Regressors,
                         Stack.Boot.MS = Stack.Boot.MS,
                         StabilitySelectionForwardKNNDetImpCoef = ForwardSelection.KNN.Det.Imp.Data.Coef,
                         StabilitySelectionForwardKNNDetImpSE = ForwardSelection.KNN.Det.Imp.Data.SE,
                         StabilitySelectionForwardKNNStocImpCoef = ForwardSelection.KNN.Stoc.Imp.Data.Coef,
                         StabilitySelectionForwardKNNStocImpSE = ForwardSelection.KNN.Stoc.Imp.Data.SE,
                         StabilitySelectionForwardRegDetImpCoef = ForwardSelection.Reg.Det.Imp.Data.Coef,
                         StabilitySelectionForwardRegDetImpSE = ForwardSelection.Reg.Det.Imp.Data.SE,
                         StabilitySelectionForwardRegStocImpCoef = ForwardSelection.Reg.Stoc.Imp.Data.Coef,
                         StabilitySelectionForwardRegStocImpSE = ForwardSelection.Reg.Stoc.Imp.Data.SE,
                         StabilitySelectionForwardMeanDetImpCoef = ForwardSelection.Mean.Det.Imp.Data.Coef,
                         StabilitySelectionForwardMeanDetImpSE = ForwardSelection.Mean.Det.Imp.Data.SE,
                         StabilitySelectionForwardMeanStocImpCoef = ForwardSelection.Mean.Stoc.Imp.Data.Coef,
                         StabilitySelectionForwardMeanStocImpSE = ForwardSelection.Mean.Stoc.Imp.Data.SE,
                         StabilitySelectionBestSubsetKNNDetImpCoef = BestSubset.KNN.Det.Imp.Data.Coef,
                         StabilitySelectionBestSubsetKNNDetImpSE = BestSubset.KNN.Det.Imp.Data.SE,
                         StabilitySelectionBestSubsetKNNStocImpCoef = BestSubset.KNN.Stoc.Imp.Data.Coef,
                         StabilitySelectionBestSubsetKNNStocImpSE = BestSubset.KNN.Stoc.Imp.Data.SE,
                         StabilitySelectionBestSubsetRegDetImpCoef = BestSubset.Reg.Det.Imp.Data.Coef,
                         StabilitySelectionBestSubsetRegDetImpSE = BestSubset.Reg.Det.Imp.Data.SE,
                         StabilitySelectionBestSubsetRegStocImpCoef = BestSubset.Reg.Stoc.Imp.Data.Coef,
                         StabilitySelectionBestSubsetRegStocImpSE = BestSubset.Reg.Stoc.Imp.Data.SE,
                         StabilitySelectionBestSubsetMeanDetImpCoef = BestSubset.Mean.Det.Imp.Data.Coef,
                         StabilitySelectionBestSubsetMeanDetImpSE = BestSubset.Mean.Det.Imp.Data.SE,
                         StabilitySelectionBestSubsetMeanStocImpCoef = BestSubset.Mean.Stoc.Imp.Data.Coef,
                         StabilitySelectionBestSubsetMeanStocImpSE = BestSubset.Mean.Stoc.Imp.Data.SE)
      )
      Sims[[MissMechanism]] <- SimResults
    }
    MissingPropResults[[MissingProportion]] <- Sims
  }
  SimulationResults[[correlation]] <- MissingPropResults
}

# The whole shebang
SimulationResults

# Save simulation results in a .Rdata file
# save(SimulationResults, file = "./Missing Data Project/MissSimulationResults.RData")

# Save all objects in the current environment to a file
# save.image("./Missing Data Project/LOO_Methods_Run.RData")