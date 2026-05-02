### Mixture Distribution Simulations
# Source functions and packages
source("./WNAR Paper/R Code/WNAR - Functions and Packages.R")

### Simulation Settings
# Seed set for each iteration within the model selection functions
nsim <- 30000

# Want plotting? If so which coefficient?
Plotting <- TRUE
PlotBeta <- "Beta5"

# Want tables?
Tabulating <- TRUE

# Sample sizes (Low / High)
n1 <- 50
n2 <- 200

# Error noise (Low / High)
eps.var1 <- 4
eps.var2 <- 16

# Information Criteria (AIC / BIC)
InfoCriteria <- c("AIC", "BIC")

# Results object
Simulation.Results <- list()
Tabulated.Results <- list()
CISubPlotDataList <- list()

index <- 0
index.ci <- 0

for(BetaSetting in c("PatternB", "PatternA")){
  # Parameter Patterns
  if (BetaSetting == "PatternA") {
    Beta <- c(0, c(0.125, 0.25, 0.50, 1, 2))
  } else {
    Beta <- c(0, c(-1.5, -0.75, 0, 0.75, 1.5))
  }
  
  # Mean vector
  mean.vec <- rep(0, length(Beta) - 1)
  
  for(SampleSize in c("LowSS", "HighSS")){
    # Establish sample size
    if (SampleSize == "LowSS") {
      sampsize <- n1
    } else {
      sampsize <- n2
    }
    
    # Establish noise level
    for(Error in c("HighError", "LowError")){
      if (Error == "LowError"){
        epsilon <- eps.var1
      } else {
        epsilon <- eps.var2
      }
      
      # Establish correlation structure
      for(Correlation in c("ModerateCorrelation", "Independent")){
        if (Correlation == "Independent") {
          Sigma1 <- diag(length(Beta) - 1)
          Sigma <- Sigma1
        } else {
          Sigma2 <- matrix(2/3, nrow = length(Beta) - 1, ncol = length(Beta) - 1)
          diag(Sigma2) <- 1
          Sigma <- Sigma2
        }
        
        # Fitting all possible models for the current simulation setting
        Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["AllModels"]] <- All.Models.Parallel.Fun(
          mu.vec = mean.vec, cov.mat = Sigma, n = sampsize, epsilon.error = epsilon, beta.seq = Beta, nsim = nsim, iterseed = TRUE)
        
        # LASSO and Ridge Model Selection
        Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][["LASSO.Ridge"]] <- Penalized.MS.Parallel.Fun(
          mu.vec = mean.vec, cov.mat = Sigma, n = sampsize, epsilon.error = epsilon,
          beta.seq = Beta, nsim = nsim, iterseed = TRUE, verbose = FALSE)
        
        if(Plotting){
          for(i in 1:(length(Beta)-1)){
            PlotBeta <- paste0("Beta", i)
          # Plot chosen coefficient (LASSO)
          Plot.Beta.Fun(Simulation.Results, MS.Type = "LASSO", SampleSize = SampleSize, Error = Error,
                        Correlation = Correlation, BetaSetting = BetaSetting, plot.stat = PlotBeta)
          
          # Plot chosen coefficient (Ridge)
          Plot.Beta.Fun(Simulation.Results, MS.Type = "Ridge", SampleSize = SampleSize, Error = Error,
                        Correlation = Correlation, BetaSetting = BetaSetting, plot.stat = PlotBeta)
          }
          stop()
        }
        if(Tabulating){
          index <- index + 1
          Tabulated.Results[[index]] <- Model.Selection.Breakdown.Fun(Simulation.Results, MS.Type = "LASSO", 
                                                                      SampleSize = SampleSize, Error = Error,
                                                                      Correlation = Correlation, IC = IC, BetaSetting = BetaSetting)
          
          index <- index + 1
          Tabulated.Results[[index]] <- Model.Selection.Breakdown.Fun(Simulation.Results, MS.Type = "Ridge", 
                                                                      SampleSize = SampleSize, Error = Error,
                                                                      Correlation = Correlation, IC = IC, BetaSetting = BetaSetting)
        }
        
        # Establish information criterion
        for(IC in InfoCriteria){
          ### Model Selection
          # Sequential Model Selection (Forward selection)
          Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][["Forward"]] <- ModelSelection.Parallel.Fun(
            mu.vec = mean.vec, cov.mat = Sigma, n = sampsize, epsilon.error = epsilon,
            beta.seq = Beta, MS.type = "forward", nsim = nsim, criteria = IC, iterseed = TRUE, verbose = FALSE)

          if(Plotting){
            # Plot chosen coefficient (Forward selection)
            for(i in 1:(length(Beta)-1)){
              PlotBeta <- paste0("Beta", i)
              Plot.Beta.Fun(Simulation.Results, MS.Type = "Forward", SampleSize = SampleSize, Error = Error,
                            Correlation = Correlation, IC = IC, BetaSetting = BetaSetting, plot.stat = PlotBeta)
            }

            Plot.Beta.Mix.Fun(Simulation.Results, MS.Type = "Forward", SampleSize = SampleSize, Error = Error,
                              Correlation = Correlation, IC = IC, BetaSetting = BetaSetting, Beta.seq = Beta)
          }

          if(Tabulating){
            index <- index + 1
            Tabulated.Results[[index]] <- Model.Selection.Breakdown.Fun(Simulation.Results, MS.Type = "Forward",
                                                                        SampleSize = SampleSize, Error = Error,
                                                                        Correlation = Correlation, IC = IC, BetaSetting = BetaSetting)
          }

          # Sequential Model Selection (Backwards selection)
          Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][["Backward"]] <- ModelSelection.Parallel.Fun(
            mu.vec = mean.vec, cov.mat = Sigma, n = sampsize, epsilon.error = epsilon,
            beta.seq = Beta, MS.type = "backward", nsim = nsim, criteria = IC, iterseed = TRUE, verbose = FALSE)

          if(Plotting){
            # Plot chosen coefficient (Backwards selection)
            Plot.Beta.Fun(Simulation.Results, MS.Type = "Backward", SampleSize = SampleSize, Error = Error,
                          Correlation = Correlation, IC = IC, BetaSetting = BetaSetting, plot.stat = PlotBeta)
          }

          if(Tabulating){
            index <- index + 1
            Tabulated.Results[[index]] <- Model.Selection.Breakdown.Fun(Simulation.Results, MS.Type = "Backward",
                                                                        SampleSize = SampleSize, Error = Error,
                                                                        Correlation = Correlation, IC = IC, BetaSetting = BetaSetting)
          }


          # Sequential Model Selection (Bi-Directional selection)
          Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][["Both"]] <- ModelSelection.Parallel.Fun(
            mu.vec = mean.vec, cov.mat = Sigma, n = sampsize, epsilon.error = epsilon,
            beta.seq = Beta, MS.type = "both", nsim = nsim, criteria = IC, iterseed = TRUE, verbose = FALSE)

          if(Plotting){
            # Plot chosen coefficient (Bi-Directional selection)
            Plot.Beta.Fun(Simulation.Results, MS.Type = "Both", SampleSize = SampleSize, Error = Error,
                          Correlation = Correlation, IC = IC, BetaSetting = BetaSetting,  plot.stat = PlotBeta)
          }

          if(Tabulating){
            index <- index + 1
            Tabulated.Results[[index]] <- Model.Selection.Breakdown.Fun(Simulation.Results, MS.Type = "Both",
                                                                        SampleSize = SampleSize, Error = Error,
                                                                        Correlation = Correlation, IC = IC, BetaSetting = BetaSetting)
          }
        }
        
        rslt.condCov <- c()
        rslt.uncondCov <- c()
        for(IC in c("AIC", "BIC")){
          for(MS.Type in c("Forward", "Backward", "Both")){
            gg <- Simulation.Results[[BetaSetting]][[SampleSize]][[Error]][[Correlation]][[IC]][[MS.Type]][["Post.Model.Selection.CIs"]]
            rslt <- CI.Coverage.Fun(gg, Beta)
            rslt.condCov <- c(rslt.condCov, rslt$Coverage.Mat[1,])
            uncondCov.temp1 <- (rslt$Coverage.Mat[1,]*rslt$Coverage.Mat[2,] + 
                                  0*(rslt$Coverage.Mat[2,1] - rslt$Coverage.Mat[2,]))/
              rslt$Coverage.Mat[2,1]
            uncondCov.temp2 <- (rslt$Coverage.Mat[1,]*rslt$Coverage.Mat[2,] + 
                                  1*(rslt$Coverage.Mat[2,1] - rslt$Coverage.Mat[2,]))/
              rslt$Coverage.Mat[2,1]
            
            uncondCov.temp <- uncondCov.temp1*(Beta != 0) + uncondCov.temp2*(Beta == 0)
            uncondCov.temp[Beta != 0 & is.na(rslt$Coverage.Mat[1,])] <- 0
            uncondCov.temp[Beta == 0 & is.na(rslt$Coverage.Mat[1,])] <- 1
            
            rslt.uncondCov <- c(rslt.uncondCov, uncondCov.temp)
          }
        }
        index.ci <- index.ci + 1
        print(paste("index.ci = ", index.ci, sep=""))
        CISubPlotDataList[[index.ci]] <- data.frame(SampleSize = SampleSize,
                                                 Error = Error,
                                                 Correlation = Correlation,
                                                 IC = rep(rep(InfoCriteria, each=18), 2),
                                                 MS.Type = rep(rep(c("Forward", "Backward", "Both"), each=6), 4),
                                                 BetaSetting = BetaSetting,
                                                 Averaging = rep(c("NAsRemoved", "NAs0"), each=36),
                                                 Beta=c("Beta0", "Beta1", "Beta2", "Beta3", "Beta4", "Beta5"),
                                                 Coverage = c(rslt.condCov, rslt.uncondCov))
        rm(Simulation.Results)
        gc()
        Simulation.Results <- list()
      }
    }
  }
}

# Save Results
# save(Tabulated.Results, file = "./WNAR Paper/TabulatedResults.RData")
# save(CISubPlotDataList, file = "./WNAR Paper/CISubPlotData.RData")

