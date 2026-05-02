### Testing Empirical Likelihood of a Univariate Mean
# Sample size, parameter dimension, 
n <- 50
d <- 1
samp <- rgamma(n, shape = 2, scale = 2)
mu0 <- 4

# Solution of weights to maximize sum of log(sample size * weights), the EL ratio,
# via method of Lagrange multipliers
wtsFun <- function(samp, lam, mu0) {
  n <- length(samp)
  wts <- (1 / n) * 1 / (1 + lam * (samp - mu0))
  return(wts)
}

# Grid search for lambda
lamSearch <- function(samp, mu0, numLam = 10000) {
  # Setting up grid given initial bounds
  n <- length(samp)
  lamMin <- (1 - (1 / n)) / (mu0 - max(samp))
  lamMax <- (1 - (1 / n)) / (mu0 - min(samp))
  lamGrid <- seq(lamMin, lamMax, length = numLam)
  
  # Calculating log empirical likelihood ratio for each lambda
  logelr <- rep(0, numLam)
  for (i in 1:numLam) {
    logelr[i] <- -2*sum(log(n * wtsFun(samp, lamGrid[i], mu0)))
  }
  
  # Return lambda EL ratio pairs
  return(data.frame(lamGrid = lamGrid, logelr = logelr))
}
lamFind <- lamSearch(samp, mu0)

# # Plotting
# plot(lamFind$lamGrid, lamFind$logelr, type = "l", xlab = "Lambda", ylab = "Log EL Ratio")
# points(lamFind[which.min(lamFind$logelr), ], col = 2, pch = 19)
# 
# # Lambda that minimizes the log empirical likelihood ratio
# optLam <- lamFind$lamGrid[which.min(lamFind$logelr)]
# 
# # Plotting weights against sample for chosen lambda
# plot(samp, wtsFun(samp, optLam, mu0), type = "h")
# # Sample mean
# abline(v = mean(samp), col = 3) 
# # Weighted mean based on chosen lambda
# abline(v = sum(wtsFun(samp, optLam, mu0) * samp), col = 4)
# # Hypothesized mean
# abline(v = mu0, col = 5)
# # Weighted mean based on chosen lambda
# sum(wtsFun(samp, optLam, mu0) * samp)

# --- 3D Plot of Weights vs Sample vs Lambda ---
# Create a lambda grid (coarser for plotting)
lamGrid <- seq(lamFind$lamGrid[1], lamFind$lamGrid[nrow(lamFind)], length.out = 1000)

# Compute weights for each lambda in the grid
wtMat <- sapply(lamGrid, function(lam) wtsFun(samp, lam, mu0))  # n x length(lamGrid)

# Create mesh for plotting
library(plotly)
x <- rep(samp, times = length(lamGrid))            # Sample value
y <- rep(lamGrid, each = length(samp))             # Lambda
z <- as.vector(wtMat)                              # Weights

# 3D scatter plot
plot_ly(
  x = x, y = y, z = z,
  type = "scatter3d", mode = "markers",
  marker = list(size = 3, color = z, colorscale = "Viridis"),
  hoverinfo = "x+y+z"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Sample X_i"),
      yaxis = list(title = "Lambda"),
      zaxis = list(title = "Weight w_i")
    ),
    title = "Empirical Likelihood Weights vs. Sample and Lambda"
  )
