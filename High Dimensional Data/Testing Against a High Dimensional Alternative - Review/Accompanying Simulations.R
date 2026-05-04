### Testing power in the F-test when we modify coefficient count and magnitude
nsim <- 10000
n <- 100
pVec <- c(5, 10, 20, 50, 90)
pVec <- seq(1, 95, by = 3)
beta0 <- 1
sigVar <- 1
betaMagVec <- seq(0, 0.25, by = 0.025)

pVals <- array(data = NA, dim = c(nsim, length(pVec), length(betaMagVec)))
for (j in 1:length(betaMagVec)) {
  pIdx <- 0
  for (p in pVec) {
    pIdx <- pIdx + 1
    betaVec <- c(beta0, rep(betaMagVec[j], p))
    for (i in 1:nsim) {
      # Data matrices
      X <- cbind(1, matrix(rnorm(n * p), nrow = n, ncol = p))
      H <- X %*% solve(t(X) %*% X) %*% t(X)
      Y <- X %*% betaVec + rnorm(n, mean = 0, sd = sqrt(sigVar))
      Yhat <- H %*% Y
      res <- Y - Yhat
      
      # Sums of squares
      SST <- (n - 1) * var(Y)
      SSE <- t(res) %*% res
      SSR <- SST - SSE
      
      # F-test
      Fstat <- (SSR / p) / (SSE / (n - p - 1))
      pVals[i, pIdx, j] <- 1 - pf(Fstat, df1 = p, df2 = n - p - 1)
    }
  }
  # Progress print
  cat("Simulations, iteration", j, "\n")
}

# Size / Power
powerMat <- colMeans(pVals < 0.05)
rownames(powerMat) <- paste0("p = ", pVec)
colnames(powerMat) <- paste0("Beta = ", betaMagVec)
powerMat

# Plotting
colvec <- rainbow(length(betaMagVec))
plot(pVec, powerMat[, 1],
     col = colvec[1],
     lwd = 1.5,
     pch = 19,
     type = "l",
     ylab = "Power",
     xlab = "p",
     main = "Power Curve \n Varying With Number Of Tested Parameters and Magnitudes", 
     ylim = c(0, 1),
     xlim = c(0, 110))
for (i in 2:length(betaMagVec)) {
  lines(pVec, powerMat[, i], col = colvec[i], lwd = 1.5)
}
legend(
  x = 97,
  y = 0.388,
  col = rev(colvec),
  lty = 1,
  lwd = 1.5,
  legend = paste0("Beta = ", rev(betaMagVec))
)

