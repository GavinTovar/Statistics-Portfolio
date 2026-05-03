library(tidyverse)
library(mvtnorm)
library(lme4)
library(gridExtra)

### Data Generation Function
Mixed.Model.DataGen <- function(
    FixEff = 4,                         # Number of fixed effects (including intercept)
    RanEff = 2,                         # Number of random effects
    Levels = rep(8, RanEff),            # Number of levels for each random effect (vector)
    N = 120,                            # Total sample size 
    Sigma.int = rep(3, RanEff),         # Std dev for random intercepts
    Sigma.slope = rep(1, RanEff),       # Std dev for random slopes 
    Cov.mat = NULL,                     # Covariance matrix for fixed-effect covariates
    Beta = NULL,                        # Fixed effect coefficients (vector)
    mu.vec = NULL,                      # Mean vector for covariates
    Response.distribution = "Normal",   # Error distribution type: "Normal", "Poisson", "Bernoulli", etc.
    Error.params = list(sd = 1),        # Parameters for error distribution (e.g., sd for normal, updated to 1)
    Link.function = "identity",         # Link function for GLMM (e.g., "identity", "log", "logit")
    overdispersion = FALSE,             # Attempt to induce overdispersion through hierarchical methods
    hierchical.str = "Poisson",         # Hierarchical structure to induce overdispersion
    Hier.params = list(lambda = 2),     # Hierarchical distribution parameters
    Noise = rnorm(N, mean = 0, sd = 2)  # Amount of noise one can such not to have a close-to-perfect-fit
) {
  # Input checks
  RanEff <- length(Levels)
  if (length(Sigma.int) != RanEff || length(Sigma.slope) != RanEff) {
    stop("Length of Sigma.int and Sigma.slope must match the number of random effects (RanEff).")
  }
  
  # Calculate sample size per level by dividing total sample size N by the number of levels in each factor
  Level.SS <- sapply(Levels, function(levels) N / levels)
  if (any(Level.SS != round(Level.SS))) {
    stop("The total sample size (N) cannot be equally divided by the number of levels. Adjust N or Levels.")
  }
  
  # Default settings
  num.cov <- FixEff - 1  # Exclude intercept
  
  # Mean vector
  if (is.null(mu.vec)) {
    mu.vec <- rep(0, num.cov)
  }
  if (length(mu.vec) != num.cov) {
    stop("mu.vec must have length equal to the number of fixed-effect covariates (FixEff - 1).")
  }
  
  # Covariance matrix
  if (is.null(Cov.mat)) {
    Cov.mat <- matrix(0.5, nrow = length(mu.vec), ncol = length(mu.vec))
    diag(Cov.mat) <- 1
  }
  if (!all(dim(Cov.mat) == num.cov)) {
    stop("Cov.mat dimensions must match the number of fixed-effect covariates (FixEff - 1).")
  }
  
  # Fixed effect beta sequence
  if (is.null(Beta)) {
    Beta <- c(1, rep(1, num.cov))
  }
  if (length(Beta) != FixEff) {
    stop("Length of Beta must equal the number of fixed effects (FixEff).")
  }
  
  # Generate fixed-effect covariates
  x.mat <- matrix(rnorm(N * num.cov, mean = mu.vec), nrow = N, ncol = num.cov)
  colnames(x.mat) <- paste0("X", 1:num.cov)
  
  # Generate random effect levels and factors
  Factors <- lapply(seq_len(RanEff), function(i) {
    rep(factor(rep(seq_len(Levels[i]), each = Level.SS[i])), length.out = N)
  })
  names(Factors) <- paste0("Z", 1:RanEff)
  
  Factors <- data.frame(Z1 = sample(Factors$Z1, N, replace = FALSE), 
                        Z2 = sample(Factors$Z2, N, replace = FALSE))
  
  # Generate random intercepts and slopes
  Random.Eff <- vector(mode = "list", length = length(RanEff))
  for(i in seq_len(RanEff)){
    intercepts <- matrix(rnorm(Levels[i], mean = 0, sd = Sigma.int[i]), nrow = Levels[i])
    slopes <- matrix(rnorm(Levels[i]*num.cov, mean = 0, sd = Sigma.slope[i]), nrow = Levels[i])
    Random.Eff[[i]] <- list(intercepts = intercepts, slopes = slopes)
  }
  # Construct the linear predictor (fixed effects + random effects)
  lin.pred <- Beta[1] + x.mat %*% Beta[-1]  # Fixed effects part
  
  # Loop through each random effect
  for (i in seq_len(RanEff)) {
    # Apply random slopes for each covariate
    for (j in seq_len(num.cov)) {
      lin.pred <- lin.pred + Random.Eff[[i]]$slopes[Factors[,i],j] * x.mat[, j]
    }
    
    # Add the random intercepts
    lin.pred <- lin.pred + Random.Eff[[i]]$intercepts[Factors[,i]]
  }
  # Apply link function (if GLMM)
  if (Link.function == "log") {
    lin.pred <- exp(lin.pred)
  } else if (Link.function == "logit") {
    lin.pred <- 1 / (1 + exp(lin.pred))
    if(any(lin.pred > 0.95)){
      print("Woah there bucko -- Inputs are creating very large success probabilities, try making them smaller if this is not intended")
    }
  }
  
  # Construct the response variable
  y.vec <- switch(Response.distribution,
                  Normal = rnorm(N, mean = lin.pred, sd = Error.params$sd),
                  Poisson = rpois(N, lambda = lin.pred),
                  Bernoulli = rbinom(N, size = 1, prob = lin.pred),
                  stop("Unsupported error distribution")
  )
  
  # Poisson overdispersion attempt
  # Thinking of the current Y values as groups and then each group has counts distributed...
  if(Response.distribution == "Poisson" & overdispersion == TRUE){
    z.vec <- rep(0, N)
    for(i in 1:N){
      if (hierchical.str == "Poisson") {
        if(!is.null(Noise)){
          z.vec[i] <- sum(rpois(y.vec[i], lambda = Hier.params$lambda)) 
          + round((Noise + Hier.params$lambda) * sqrt(Hier.params$lambda))
        }
        z.vec[i] <- sum(rpois(y.vec[i], lambda = Hier.params$lambda))
      } else if (hierchical.str == "Gamma") {
        if(!is.null(Noise)){
          z.vec[i] <- round(sum(rgamma(y.vec[i], shape = Hier.params$shape, rate = Hier.params$rate))  
                            + (Noise + Hier.params$shape/Hier.params$shape) * sqrt(Hier.params$shape/(Hier.params$shape)^2))
        }
        z.vec[i] <- round(sum(rgamma(y.vec[i], shape = Hier.params$shape, rate = Hier.params$rate))) 
      } else {
        if(!is.null(Noise)){
          z.vec[i] <- round(sum(rnorm(y.vec[i], mean = Hier.params$mean, sd = Hier.params$sd))
                            + (Noise + Hier.params$mean) * Hier.params$sd)
        }
        z.vec[i] <- round(sum(rnorm(y.vec[i], mean = Hier.params$mean, sd = Hier.params$sd)))
      }
    }
    y.vec <- z.vec
  }
  
  if(!is.null(Noise) & Response.distribution != "Bernoulli" & overdispersion == FALSE){
    if(Response.distribution == "Poisson") {
      y.vec <- y.vec + round(abs(Noise))
    } else {
      y.vec <- y.vec + Noise
    }
    
  }
  
  # Combine data into a dataframe
  df <- data.frame(Y = y.vec, x.mat, do.call(cbind, Factors))
  return(df)
}

# Simulations
# Standard GLMM data (no overdispersion)
set.seed(10262025)
nsim <- 300

# Object Creation
AIC.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Dev.Mat <- matrix(data = NA, nrow = nsim, ncol = 4)
Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Ran1Eff <- Ran2Eff <- matrix(data = NA, nrow = nsim, ncol = 3)

for(i in 1:nsim){
  # Generate Data
  sim.data0 <- Mixed.Model.DataGen(
    FixEff = 3,
    RanEff = 2,
    Levels = c(40, 60),
    N = 2400,
    Sigma.int = c(0.50, 0.50),
    Sigma.slope = c(1, 1),
    Cov.mat = NULL,
    Beta = c(0.25, 0.50, 0.75),
    mu.vec = c(0, 0),
    Response.distribution = "Poisson",
    Error.params = NULL,
    Link.function = "log",
    overdispersion = FALSE,
    Noise = NULL
  )
  
  # Attempt to fit the model
  PoisMod0 <- tryCatch({
    glmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
          data = sim.data0, 
          family = poisson(link = "log"))
  }, error = function(e) {
    message("Model failed to fit: ", e$message)
    return(NULL) # Return NULL if an error occurs
  })
  
  # Check if the model fitting failed
  if (is.null(PoisMod0)) {
    next
  }
  
  # Model Summary
  sum.mod <- summary(PoisMod0)
  
  # Pull AIC/BIC/Log Likelihood
  AIC.Mat[i, ] <- sum.mod$AICtab[1:3]
  
  # Residual Deviance, Deviance, DF Residual, OD Stat
  Dev.Mat[i, ] <- c(sum.mod$AICtab[4], deviance(PoisMod0), 
                    sum.mod$AICtab[5], deviance(PoisMod0)/sum.mod$AICtab[5])
  
  # Extracting Fixed Effect Statistics
  Est.Mat[i, ] <- sum.mod$coefficients[, 1]
  SE.Est.Mat[i, ] <- sum.mod$coefficients[, 2]
  Z.Stat.Mat[i, ] <- sum.mod$coefficients[, 3]
  
  # Extracting Random Effects
  Ran1Eff[i, ] <- attr(sum.mod$varcor[[1]], "stddev")
  Ran2Eff[i, ] <- attr(sum.mod$varcor[[2]], "stddev")
  
  # Progress Print
  if(i %% 10 == 0) print(i)
  
}
# Labeling
colnames(AIC.Mat) <- c("AIC", "BIC", "LogLik")
colnames(Dev.Mat) <- c("Resid.Dev.", "Dev.", "df.Resid", "OD")
colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2")
colnames(Ran1Eff) <- paste0("Z1", c("(Intercept)", "X1", "X2"))
colnames(Ran2Eff) <- paste0("Z2", c("(Intercept)", "X1", "X2"))


### Simple Plotting Function
Plot.Thing <- function(data, plot.item, Color, True.Coef = NA) {
  plot <- ggplot(data, aes(x = .data[[plot.item]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = Color, col = "black", alpha = 0.25) +
    geom_density(color = Color, linewidth = 1.2, alpha = 0.66) +
    labs(
      x = plot.item, 
      y = "Density", 
      title = paste("Density Plot of", plot.item),
      subtitle = "Histogram Overlaid with Density Curve"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(face = "italic", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if(!is.na(True.Coef)){
    plot <- plot + geom_vline(xintercept = True.Coef, color = "black", linetype = "dashed", linewidth = 1)
  }
  return(plot)
}

# Plotting Random Effect Estimates
p1 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1(Intercept)", Color = "#68fc8f", True.Coef = 0.50)
p2 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X1", Color = "#42fc73", True.Coef = 1)
p3 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X2", Color = "#17ff55", True.Coef = 1)
p4 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2(Intercept)", Color = "#83b4fc", True.Coef = 0.50)
p5 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X1", Color = "#599cff", True.Coef = 1)
p6 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X2", Color = "#2179fc", True.Coef = 1)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# Plotting Fixed Effect Estimates
p1 <- Plot.Thing(data = Est.Mat, plot.item = "(Intercept)", Color = "#fc77f8", True.Coef = 0.25)
p2 <- Plot.Thing(data = Est.Mat, plot.item = "X1", Color = "#fc4cf7", True.Coef = 0.50)
p3 <- Plot.Thing(data = Est.Mat, plot.item = "X2", Color = "#ff26f8", True.Coef = 0.75)
grid.arrange(p1, p2, p3, nrow = 1)

### Getting some stats
apply(AIC.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Dev.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran1Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran2Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)


### Some takeaways from the above
# Settings 1
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(4, 5),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Settings 2 (Increase SS)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(4, 5),
  N = 2400,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Settings 3 (Increase level counts)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 50),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Settings 4 (Increase SS and level counts)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 50),
  N = 2400,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Low random effect coefficients increase the probability of singularity
# Very small random effects can result in convergence problems
# When random effect levels counts are low we see low accuracy and precision in the fixed effect estimates
# 10x sample size doesn't meaningfully decrease fixed effect SE or RE estimates

### Exploration 2
# Fitting a GLM to GLMM data - So ignoring the random effects present
set.seed(10262025)
nsim <- 50000

# Object Creation
Info.Mat <- matrix(data = NA, nrow = nsim, ncol = 4)
Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 101)
# Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 24)
# Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 9)

for(i in 1:nsim){
  # Generate Data
  sim.data0 <- Mixed.Model.DataGen(
    FixEff = 3,
    RanEff = 2,
    Levels = c(40, 60),
    N = 240,
    Sigma.int = c(0.50, 0.50),
    Sigma.slope = c(1, 1),
    Cov.mat = NULL,
    Beta = c(0.25, 0.50, 0.75),
    mu.vec = c(0, 0),
    Response.distribution = "Poisson",
    Error.params = NULL,
    Link.function = "log",
    overdispersion = FALSE,
    Noise = NULL
  )
  
  # Attempt to fit the model with error and warning handling
  PoisGLM <- tryCatch({
    withCallingHandlers(
      glm(Y ~ X1 + X2
            # + X1*Z1 + X1*Z2 + X2*Z1 + X2*Z2,
            + factor(Z1) + factor(Z2),
            #factor(Z1)*X1 + factor(Z1)*X2 + factor(Z2)*X1 + factor(Z2)*X2, 
          data = sim.data0,
          family = poisson(link = "log")),
      warning = function(w) {
        # Check for specific warnings
        if (grepl("fitted rates numerically 0 occurred|algorithm did not converge", w$message)) {
          message("Warning caught: ", w$message)
          invokeRestart("muffleWarning")
          stop("Skipping due to warning: ", w$message)
        }
      }
    )
  }, error = function(e) {
    message("Model failed to fit: ", e$message)
    return(NULL)
  })
  
  # Skip iteration if the model fitting failed
  if (is.null(PoisGLM)) {
    next
  }
  
  # Model Summary
  sum.mod <- summary(PoisGLM)
  
  # Residual Deviance, DF Residual, OD Stat, AIC
  Info.Mat[i, ] <- c(sum.mod$deviance, sum.mod$df.residual, 
                     sum.mod$deviance/sum.mod$df.residual, sum.mod$aic)
  
  # Skip if any coeffs are NA
  if(any(sum.mod$aliased)){
    next
  }
  
  # Extracting Fixed Effect Statistics
  Est.Mat[i, ] <- sum.mod$coefficients[, 1]
  SE.Est.Mat[i, ] <- sum.mod$coefficients[, 2]
  Z.Stat.Mat[i, ] <- sum.mod$coefficients[, 3]
  
  # Progress Print
  if(i %% 1000 == 0) print(i)
  
}

# Labeling
colnames(Info.Mat) <- c("Resid.Dev.", "df.Resid", "OD", "AIC")
# colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2", "Z1", "Z2",
#                                                                        "X1:Z1", "X1:Z2", "X2:Z1", "X2:Z2")
# colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2", 
#                                                                        paste0("(Z1)", 2:4), paste0("(Z2)", 2:5))
colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2",
                                                                       paste0("(Z1)", 2:40), paste0("(Z2)", 2:60))
# colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2", 
#                                                                        paste0("(Z1)", 2:4), paste0("(Z2)", 2:5),
#                                                                        paste0("(X1:Z1)", 2:4), paste0("(X2:Z1)", 2:4),
#                                                                        paste0("(X1:Z2)", 2:5), paste0("(X2:Z2)", 2:5))

### Getting some stats
apply(Info.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)

# Table for LaTeX
rbind.data.frame(t(apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)),
                 t(apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)),
                 t(apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)))

result <- rbind.data.frame(
  t(apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)),
  t(apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)),
  t(apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE))
)
result[] <- lapply(result, function(x) format(x, scientific = FALSE))
print(result)


### Simple Plotting Function
Plot.Thing <- function(data, plot.item, Color, True.Coef = NA) {
  plot <- ggplot(data, aes(x = .data[[plot.item]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = Color, col = "black", alpha = 0.25) +
    geom_density(color = Color, linewidth = 1.2, alpha = 0.66) +
    labs(
      x = plot.item, 
      y = "Density", 
      title = paste("Density Plot of", plot.item),
      subtitle = "Histogram Overlaid with Density Curve"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(face = "italic", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if(!is.na(True.Coef)){
    plot <- plot + geom_vline(xintercept = True.Coef, color = "black", linetype = "dashed", linewidth = 1)
  }
  return(plot)
}

# Plotting
# Intercept - X1 - X2
p0 <- Plot.Thing(data = Est.Mat, plot.item = "(Intercept)", Color = "#A8E6A3", True.Coef = 0.25)
p1 <- Plot.Thing(data = Est.Mat, plot.item = "X1", Color = "#66C66C", True.Coef = 0.50)
p2 <- Plot.Thing(data = Est.Mat, plot.item = "X2", Color = "#006400", True.Coef = 0.75)
grid.arrange(p0, p1, p2, ncol = 3)

# Random Effects
# Define the inputs
plot_items1 <- c(paste0("(Z1)", 2:4), paste0("(Z2)", 2:5))

# Define a range of pink shades
pink_shades1 <- c("#FFC6DD", "#FFB6D5", "#FFA7CE", "#FF97C6", "#FF88BE", 
                 "#FF78B7", "#FF69AF")

# Create plots in a loop, assigning each plot a different shade of pink
plots <- list()

for (i in 1:length(plot_items1)) {
  # Dynamically create the plot object using Plot.Thing
  plot_name <- paste0("p", i)
  plots[[plot_name]] <- Plot.Thing(data = Est.Mat, 
                                   plot.item = plot_items1[i], 
                                   Color = pink_shades1[i], 
                                   True.Coef = 0)
}
layout1 <- rbind(c(1, 1, 2, 2, 3, 3),
                 c(4, 4, 4, 5, 5, 5),
                 c(6, 6, 6, 7, 7, 7))
grid.arrange(grobs = plots, ncol = 3, layout_matrix = layout1)

# Random Effects and interactions
# Define the inputs
plot_items2 <- c(paste0("(X1:Z1)", 2:4), paste0("(X2:Z1)", 2:4),
                 paste0("(X1:Z2)", 2:5), paste0("(X2:Z2)", 2:5))

# Define a range of pink shades
pink_shades2 <- c("#FF59A7", "#FF4A9F", "#FF3A97", 
                  "#FF2B8F", "#FF1B87", "#FF0C7F", "#F50077", "#E0006F", 
                  "#CC0067", "#B8005F", "#A30057", "#8F004F", "#7A0047", "#66003F")

# Create plots in a loop, assigning each plot a different shade of pink
plots2 <- list()

for (i in 1:length(plot_items2)) {
  # Dynamically create the plot object using Plot.Thing
  plot_name <- paste0("p", i)
  plots2[[plot_name]] <- Plot.Thing(data = Est.Mat, 
                                   plot.item = plot_items2[i], 
                                   Color = pink_shades2[i], 
                                   True.Coef = 0)
}

layout1 <- rbind(rep(1:3, each = 4),
                 rep(4:6, each = 4),
                 rep(7:10, each = 3),
                 rep(11:14, each = 3))

grid.arrange(grobs = plots2, nrow = 4, layout_matrix = layout1)


# Settings 1
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(4, 5),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Settings 2 (Increase SS)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(4, 5),
  N = 2400,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Settings 3 (Increase level counts)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 50),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)

# Cannot fit model with this many levels. Coefficients: (57 not defined because of singularities)

### Exploration 3
# Exploring overdispersion through hierarchical means
set.seed(10262025)
nsim <- 300

# Object Creation
AIC.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Dev.Mat <- matrix(data = NA, nrow = nsim, ncol = 4)
Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Ran1Eff <- Ran2Eff <- matrix(data = NA, nrow = nsim, ncol = 3)

for(i in 1:nsim){
  # Generate Data
  sim.data0 <- Mixed.Model.DataGen(
    FixEff = 3,
    RanEff = 2,
    Levels = c(40, 60),
    N = 240,
    Sigma.int = c(0.50, 0.50),
    Sigma.slope = c(1, 1),
    Cov.mat = NULL,
    Beta = c(0.25, 0.50, 0.75),
    mu.vec = c(0, 0),
    Response.distribution = "Poisson",
    Error.params = NULL,
    Link.function = "log",
    overdispersion = TRUE,
    hierchical.str = "Normal",
    Hier.params = list(mean = 6, sd = 2),
    Noise = NULL
  )
  
  # Attempt to fit the model
  PoisMod0 <- tryCatch({
    glmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
          data = sim.data0, 
          family = poisson(link = "log"))
  }, error = function(e) {
    message("Model failed to fit: ", e$message)
    return(NULL) # Return NULL if an error occurs
  })
  
  # Check if the model fitting failed
  if (is.null(PoisMod0)) {
    next
  }
  
  # Model Summary
  sum.mod <- summary(PoisMod0)
  
  # Pull AIC/BIC/Log Likelihood
  AIC.Mat[i, ] <- sum.mod$AICtab[1:3]
  
  # Residual Deviance, Deviance, DF Residual, OD Stat
  Dev.Mat[i, ] <- c(sum.mod$AICtab[4], deviance(PoisMod0), 
                    sum.mod$AICtab[5], deviance(PoisMod0)/sum.mod$AICtab[5])
  
  # Extracting Fixed Effect Statistics
  Est.Mat[i, ] <- sum.mod$coefficients[, 1]
  SE.Est.Mat[i, ] <- sum.mod$coefficients[, 2]
  Z.Stat.Mat[i, ] <- sum.mod$coefficients[, 3]
  
  # Extracting Random Effects
  Ran1Eff[i, ] <- attr(sum.mod$varcor[[1]], "stddev")
  Ran2Eff[i, ] <- attr(sum.mod$varcor[[2]], "stddev")
  
  # Progress Print
  if(i %% 10 == 0) print(i)
  
}
# Labeling
colnames(AIC.Mat) <- c("AIC", "BIC", "LogLik")
colnames(Dev.Mat) <- c("Resid.Dev.", "Dev.", "df.Resid", "OD")
colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2")
colnames(Ran1Eff) <- paste0("Z1", c("(Intercept)", "X1", "X2"))
colnames(Ran2Eff) <- paste0("Z2", c("(Intercept)", "X1", "X2"))


### Simple Plotting Function
Plot.Thing <- function(data, plot.item, Color, True.Coef = NA) {
  plot <- ggplot(data, aes(x = .data[[plot.item]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = Color, col = "black", alpha = 0.25) +
    geom_density(color = Color, linewidth = 1.2, alpha = 0.66) +
    labs(
      x = plot.item, 
      y = "Density", 
      title = paste("Density Plot of", plot.item),
      subtitle = "Histogram Overlaid with Density Curve"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(face = "italic", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if(!is.na(True.Coef)){
    plot <- plot + geom_vline(xintercept = True.Coef, color = "black", linetype = "dashed", linewidth = 1)
  }
  return(plot)
}

# Plotting Random Effect Estimates
p1 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1(Intercept)", Color = "#68fc8f", True.Coef = 0.50)
p2 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X1", Color = "#42fc73", True.Coef = 1)
p3 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X2", Color = "#17ff55", True.Coef = 1)
p4 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2(Intercept)", Color = "#83b4fc", True.Coef = 0.50)
p5 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X1", Color = "#599cff", True.Coef = 1)
p6 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X2", Color = "#2179fc", True.Coef = 1)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# Plotting Fixed Effect Estimates
p1 <- Plot.Thing(data = Est.Mat, plot.item = "(Intercept)", Color = "#fc77f8", True.Coef = 0.25)
p2 <- Plot.Thing(data = Est.Mat, plot.item = "X1", Color = "#fc4cf7", True.Coef = 0.50)
p3 <- Plot.Thing(data = Est.Mat, plot.item = "X2", Color = "#ff26f8", True.Coef = 0.75)
grid.arrange(p1, p2, p3, nrow = 1)

### Getting some stats
apply(AIC.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Dev.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran1Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran2Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)

# Setting 1
# Same data generation but now with overdispersion via Poisson hierarchy
Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 60),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = TRUE,
  hierchical.str = "Poisson",
  Hier.params = list(lambda = 3),
  Noise = NULL
)

# Setting 2
# Same data generation but now with overdispersion via Gamma hierarchy
Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 60),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = TRUE,
  hierchical.str = "Gamma",
  Hier.params = list(shape = 5, rate = 3),
  Noise = NULL
)


# Setting 3
# Same data generation but now with overdispersion via Normal hierarchy
Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 60),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = TRUE,
  hierchical.str = "Normal",
  Hier.params = list(mean = 6, sd = 2),
  Noise = NULL
)


### Explortation 4
# Different levels of noise and then
set.seed(10262025)
nsim <- 300

# Object Creation
AIC.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Dev.Mat <- matrix(data = NA, nrow = nsim, ncol = 4)
Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
Ran1Eff <- Ran2Eff <- matrix(data = NA, nrow = nsim, ncol = 3)

for(i in 1:nsim){
  # Generate Data
  sim.data0 <- Mixed.Model.DataGen(
    FixEff = 3,
    RanEff = 2,
    Levels = c(40, 60),
    N = 2400,
    Sigma.int = c(0.50, 0.50),
    Sigma.slope = c(1, 1),
    Cov.mat = NULL,
    Beta = c(0.25, 0.50, 0.75),
    mu.vec = c(0, 0),
    Response.distribution = "Poisson",
    Error.params = NULL,
    Link.function = "log",
    overdispersion = FALSE,
    Noise = NULL
  )
  
  # Attempt to fit the model
  PoisMod0 <- tryCatch({
    glmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2),
          data = sim.data0, 
          family = poisson(link = "log"))
  }, error = function(e) {
    message("Model failed to fit: ", e$message)
    return(NULL) # Return NULL if an error occurs
  })
  
  # Check if the model fitting failed
  if (is.null(PoisMod0)) {
    next
  }
  
  # Model Summary
  sum.mod <- summary(PoisMod0)
  
  # Pull AIC/BIC/Log Likelihood
  AIC.Mat[i, ] <- sum.mod$AICtab[1:3]
  
  # Residual Deviance, Deviance, DF Residual, OD Stat
  Dev.Mat[i, ] <- c(sum.mod$AICtab[4], deviance(PoisMod0), 
                    sum.mod$AICtab[5], deviance(PoisMod0)/sum.mod$AICtab[5])
  
  # Extracting Fixed Effect Statistics
  Est.Mat[i, ] <- sum.mod$coefficients[, 1]
  SE.Est.Mat[i, ] <- sum.mod$coefficients[, 2]
  Z.Stat.Mat[i, ] <- sum.mod$coefficients[, 3]
  
  # Extracting Random Effects
  Ran1Eff[i, ] <- attr(sum.mod$varcor[[1]], "stddev")
  Ran2Eff[i, ] <- attr(sum.mod$varcor[[2]], "stddev")
  
  # Progress Print
  if(i %% 10 == 0) print(i)
  
}
# Labeling
colnames(AIC.Mat) <- c("AIC", "BIC", "LogLik")
colnames(Dev.Mat) <- c("Resid.Dev.", "Dev.", "df.Resid", "OD")
colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2")
colnames(Ran1Eff) <- paste0("Z1", c("(Intercept)", "X1", "X2"))
colnames(Ran2Eff) <- paste0("Z2", c("(Intercept)", "X1", "X2"))


### Simple Plotting Function
Plot.Thing <- function(data, plot.item, Color, True.Coef = NA) {
  plot <- ggplot(data, aes(x = .data[[plot.item]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = Color, col = "black", alpha = 0.25) +
    geom_density(color = Color, linewidth = 1.2, alpha = 0.66) +
    labs(
      x = plot.item, 
      y = "Density", 
      title = paste("Density Plot of", plot.item),
      subtitle = "Histogram Overlaid with Density Curve"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(face = "italic", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if(!is.na(True.Coef)){
    plot <- plot + geom_vline(xintercept = True.Coef, color = "black", linetype = "dashed", linewidth = 1)
  }
  return(plot)
}

# Plotting Random Effect Estimates
p1 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1(Intercept)", Color = "#68fc8f", True.Coef = 0.50)
p2 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X1", Color = "#42fc73", True.Coef = 1)
p3 <- Plot.Thing(data = Ran1Eff, plot.item = "Z1X2", Color = "#17ff55", True.Coef = 1)
p4 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2(Intercept)", Color = "#83b4fc", True.Coef = 0.50)
p5 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X1", Color = "#599cff", True.Coef = 1)
p6 <- Plot.Thing(data = Ran2Eff, plot.item = "Z2X2", Color = "#2179fc", True.Coef = 1)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# Plotting Fixed Effect Estimates
p1 <- Plot.Thing(data = Est.Mat, plot.item = "(Intercept)", Color = "#fc77f8", True.Coef = 0.25)
p2 <- Plot.Thing(data = Est.Mat, plot.item = "X1", Color = "#fc4cf7", True.Coef = 0.50)
p3 <- Plot.Thing(data = Est.Mat, plot.item = "X2", Color = "#ff26f8", True.Coef = 0.75)
grid.arrange(p1, p2, p3, nrow = 1)

### Getting some stats
apply(AIC.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Dev.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran1Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Ran2Eff, MARGIN = 2, FUN = mean, na.rm = TRUE)

# Setting 1
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 60),
  N = 240,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = rnorm(240, mean = 0, sd = 2)
)

# Setting 2
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(40, 60),
  N = 2400,
  Sigma.int = c(0.50, 0.50),
  Sigma.slope = c(1, 1),
  Cov.mat = NULL,
  Beta = c(0.25, 0.50, 0.75),
  mu.vec = c(0, 0),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = rnorm(2400, mean = 0, sd = 4)
)



### GLM w/ no random effects
set.seed(10262025)
nsim <- 50000

# Object Creation
Info.Mat <- matrix(data = NA, nrow = nsim, ncol = 4)
Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 3)
# Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 101)
# Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 24)
# Est.Mat <- SE.Est.Mat <- Z.Stat.Mat <- matrix(data = NA, nrow = nsim, ncol = 9)

for(i in 1:nsim){
  # Generate Data
  sim.data0 <- Mixed.Model.DataGen(
    FixEff = 3,
    RanEff = 2,
    Levels = c(4, 5),
    N = 240,
    Sigma.int = c(0.50, 0.50),
    Sigma.slope = c(1, 1),
    Cov.mat = NULL,
    Beta = c(0.25, 0.50, 0.75),
    mu.vec = c(0, 0),
    Response.distribution = "Poisson",
    Error.params = NULL,
    Link.function = "log",
    overdispersion = FALSE,
    Noise = NULL
  )
  
  # Attempt to fit the model with error and warning handling
  PoisGLM <- tryCatch({
    withCallingHandlers(
      glm(Y ~ X1 + X2, 
          data = sim.data0,
          family = poisson(link = "log")),
      warning = function(w) {
        # Check for specific warnings
        if (grepl("fitted rates numerically 0 occurred|algorithm did not converge", w$message)) {
          message("Warning caught: ", w$message)
          invokeRestart("muffleWarning")
          stop("Skipping due to warning: ", w$message)
        }
      }
    )
  }, error = function(e) {
    message("Model failed to fit: ", e$message)
    return(NULL)
  })
  
  # Skip iteration if the model fitting failed
  if (is.null(PoisGLM)) {
    next
  }
  
  # Model Summary
  sum.mod <- summary(PoisGLM)
  
  # Residual Deviance, DF Residual, OD Stat, AIC
  Info.Mat[i, ] <- c(sum.mod$deviance, sum.mod$df.residual, 
                     sum.mod$deviance/sum.mod$df.residual, sum.mod$aic)
  
  # Skip if any coeffs are NA
  if(any(sum.mod$aliased)){
    next
  }
  
  # Extracting Fixed Effect Statistics
  Est.Mat[i, ] <- sum.mod$coefficients[, 1]
  SE.Est.Mat[i, ] <- sum.mod$coefficients[, 2]
  Z.Stat.Mat[i, ] <- sum.mod$coefficients[, 3]
  
  # Progress Print
  if(i %% 1000 == 0) print(i)
  
}

# Labeling
colnames(Info.Mat) <- c("Resid.Dev.", "df.Resid", "OD", "AIC")
colnames(Est.Mat) <- colnames(SE.Est.Mat) <- colnames(Z.Stat.Mat) <- c("(Intercept)", "X1", "X2")


### Simple Plotting Function
Plot.Thing <- function(data, plot.item, Color, True.Coef = NA) {
  plot <- ggplot(data, aes(x = .data[[plot.item]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = Color, col = "black", alpha = 0.25) +
    geom_density(color = Color, linewidth = 1.2, alpha = 0.66) +
    labs(
      x = plot.item, 
      y = "Density", 
      title = paste("Density Plot of", plot.item),
      subtitle = "Histogram Overlaid with Density Curve"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(face = "italic", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if(!is.na(True.Coef)){
    plot <- plot + geom_vline(xintercept = True.Coef, color = "black", linetype = "dashed", linewidth = 1)
  }
  return(plot)
}

# Plotting Fixed Effect Estimates
p1 <- Plot.Thing(data = Est.Mat, plot.item = "(Intercept)", Color = "#fc77f8", True.Coef = 0.25)
p2 <- Plot.Thing(data = Est.Mat, plot.item = "X1", Color = "#fc4cf7", True.Coef = 0.50)
p3 <- Plot.Thing(data = Est.Mat, plot.item = "X2", Color = "#ff26f8", True.Coef = 0.75)
grid.arrange(p1, p2, p3, nrow = 1)

### Getting some stats
apply(Info.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(SE.Est.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)
apply(Z.Stat.Mat, MARGIN = 2, FUN = mean, na.rm = TRUE)