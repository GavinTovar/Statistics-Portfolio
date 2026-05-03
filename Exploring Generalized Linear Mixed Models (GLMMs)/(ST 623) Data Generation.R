### Set seed for reproducibility and load required packages
set.seed(10262025)
library(tidyverse)
library(mvtnorm)
library(lme4)
library(patchwork)

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
  
  # Generate random intercepts and slopes
  Random.Eff <- lapply(seq_len(RanEff), function(i) {
    intercepts <- rnorm(Levels[i], mean = 0, sd = Sigma.int[i])
    slopes <- if (Sigma.slope[i] > 0) {
      rnorm(Levels[i], mean = 0, sd = Sigma.slope[i])
      } else {
        rep(0, Levels[i])
      }
    list(intercepts = intercepts[as.integer(Factors[[i]])], slopes = slopes[as.integer(Factors[[i]])])
  })
  
  
  #####
  
  ### Perhaps unneeded randomization element???
  # Randomize order of random effect variables and align Factors
  Random.Eff <- lapply(seq_along(Random.Eff), function(i) {
    # Randomize order of intercepts and slopes
    random_indices <- sample(seq_along(Random.Eff[[i]]$intercepts))
    list(
      intercepts = Random.Eff[[i]]$intercepts[random_indices],
      slopes = Random.Eff[[i]]$slopes[random_indices],
      random_indices = random_indices # Store the order for aligning Factors
    )
  })
  
  # Update Factors to match the randomized order
  Factors <- lapply(seq_along(Factors), function(i) {
    factor(as.integer(Factors[[i]])[Random.Eff[[i]]$random_indices])
  })
  names(Factors) <- paste0("Z", seq_along(Factors))
  
  #####
  
  # Construct the linear predictor (fixed effects + random effects)
  lin.pred <- Beta[1] + x.mat %*% Beta[-1]  # Fixed effects part

  # Loop through each random effect
  for (i in seq_len(RanEff)) {
    # Apply random slopes for each covariate
    for (j in seq_len(num.cov)) {
      lin.pred <- lin.pred + Random.Eff[[i]]$slopes[j] * x.mat[, j]
    }
    
    # Add the random intercepts
    lin.pred <- lin.pred + Random.Eff[[i]]$intercepts
  }
  browser()
  # Apply link function (if GLMM)
  if (Link.function == "log") {
    lin.pred <- exp(lin.pred)
    if(any(lin.pred > 100)){
      print("Woah there bucko -- Inputs are creating very large counts, try making them smaller if this is not intended")
    }
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

# Example: Default settings - Simulating Normal data (LMM data)
sim.data <- Mixed.Model.DataGen()
sim.data

# Fit GLMM with exact model
NormMod <- glmer(Y ~ X1 + X2 + X3 + (X1 + X2 + X3 | Z1) + (X1 + X2 + X3 | Z2), data = sim.data)
summary(NormMod)

# Simulate data with adjusted parameters
sim.data1.2 <- Mixed.Model.DataGen(N = 120, 
                                   Beta = c(1, 2, 2, 2),
                                   Sigma.int = rep(3, 2), 
                                   Error.params = list(sd = 1))
# Fit GLMM
NormMod2 <- lmer(Y ~ X1 + X2 + X3 + (X1 | Z1) + (X1 | Z2), data = sim.data1.2)
summary(NormMod2)

# Example: Simulating Bernoulli GLMM data
sim.data2 <- Mixed.Model.DataGen(
  FixEff = 4,
  RanEff = 3,
  Levels = c(6, 10, 4),
  N = 120,
  Sigma.int = c(0.01, 0.05, 0.02),
  Sigma.slope = c(0.5, 0, 0.01),
  Cov.mat = matrix(c(1, 0.2, 0.2, 1, 0.5, 0.5, 0.3, 0.25, 0.4), 3, 3),
  Beta = c(0.04, 0.03, 0.04, 0.01),
  mu.vec = c(0, 1, 1),
  Response.distribution = "Bernoulli",
  Error.params = NULL,
  Link.function = "logit"
)
sim.data2

# Fit GLMM
BernMod <- glmer(Y ~ X1 + X2 + X3 + (1 | Z1) + (1 | Z2) + (1 | Z3), 
                 data = sim.data2, 
                 family = binomial(link = "logit"))

# Model summary
summary(BernMod)

# Simulating Poisson GLMM data
# Keep coefficients < 0.66 as a rule of thumb - or big counts might result
sim.data3 <- Mixed.Model.DataGen(
  FixEff = 4,
  RanEff = 3,
  Levels = c(6, 5, 4),
  N = 240,
  Sigma.int = c(0.05, 0.05, 0.05),
  Sigma.slope = c(0.15, 0.15, 0.15),
  Cov.mat = NULL,
  Beta = c(0.15, 0.25, 0.50, 0.66),
  mu.vec = c(1, 1, 1.25),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE
)
sim.data3

# Fit GLMM
PoisMod <- glmer(Y ~ X1 + X2 + X3 + (X1 + X2 + X3 | Z1) + (X1 + X2 + X3 | Z2) + (X1 + X2 + X3 | Z3), 
             data = sim.data3, 
             family = poisson(link = "log"))

# Model summary
summary(PoisMod)

#####
### Random Effect Plotting Function (A work in progress)
plot_random_effects_vs_fixed <- function(model, 
                                         data, 
                                         fixed_effects = NULL, 
                                         random_effect = NULL, 
                                         include_slopes = TRUE, 
                                         plot_response = FALSE, 
                                         facet = FALSE, 
                                         plot_summary_stats = FALSE) {
  
  # Extract random effects from the GLMM model
  ranef_data <- ranef(model, condVar = TRUE)
  
  # Extract the fixed effects and predictions from the model
  data$Predicted <- predict(model, newdata = data, re.form = NA) # Fixed effects predictions
  
  # If plotting summary statistics is requested
  if (plot_summary_stats) {
    # Initialize an empty list to store the random effects data
    ranef_df_list <- list()
    
    # Loop over each random effect group in the model
    for (group_name in names(ranef_data)) {
      # Extract random intercept and slope for the group
      random_intercept <- ranef_data[[group_name]][, "(Intercept)"]
      
      # If the group has random slopes, extract them as well
      random_slopes <- ifelse(ncol(ranef_data[[group_name]]) > 1, ranef_data[[group_name]][, 2], NA)
      
      # Prepare data for plotting
      group_data <- data.frame(
        Group = rownames(ranef_data[[group_name]]),
        RandomIntercept = random_intercept,
        RandomSlope = random_slopes,
        Factor = group_name
      )
      
      # Append to the list
      ranef_df_list[[group_name]] <- group_data
    }
    
    # Combine all random effects data into one data frame
    ranef_df <- do.call(rbind, ranef_df_list)
    
    # Check and ensure proper matching of data columns for fixed effects
    # fixed_effects_summary <- data.frame(Z1 = unique(data$Z1))
    fixed_effects_summary <- data.frame(setNames(list(unique(data[[random_effect]])), random_effect))
    
    # Add mean values for each fixed effect to the fixed_effects_summary
    for (fixed in fixed_effects) {
      fixed_mean <- aggregate(as.formula(paste(fixed, "~", random_effect)), data = data, FUN = mean)
      fixed_effects_summary <- merge(fixed_effects_summary, fixed_mean, by = random_effect, all.x = TRUE)
    }
    
    # Merge random effects with fixed effects summary
    merged_data <- merge(ranef_df, fixed_effects_summary, by.x = "Group", by.y = random_effect, all.x = TRUE)
    
    # Create the plot for each fixed effect
    plot_list <- lapply(seq_along(fixed_effects), function(index) {
      fixed_effect <- fixed_effects[index]
      
      plot <- ggplot(merged_data, aes_string(x = fixed_effect, y = "RandomIntercept", color = "Factor")) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "black") +
        geom_line(aes(y = Predicted), linetype = "dashed") +  # Add model-predicted regression line
        labs(
          title = paste("Random Effects vs", fixed_effect),
          x = fixed_effect,
          y = "Random Intercept",
          color = "Random Effects"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      # Remove the Y-axis elements for all but the first plot
      if (index > 1) {
        plot <- plot + theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
      
      return(plot)
    })
    
    # Combine plots with patchwork
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = length(fixed_effects)) + 
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    # Return the combined plot
    return(combined_plot)
    
  } else {
    # Extract specified random effect data
    if (!is.null(random_effect) && random_effect %in% names(ranef_data)) {
      random_effect_data <- ranef_data[[random_effect]]
      
      # Check if the model has random slopes
      if (include_slopes && ncol(random_effect_data) > 1) {
        intercepts <- random_effect_data[, "(Intercept)"]
        slopes <- ifelse(ncol(random_effect_data) > 1, random_effect_data[, 2], NA)  # Dynamically use the fixed effect variable
        
        # Prepare data for plotting
        ranef_df <- data.frame(
          Group = rownames(random_effect_data),
          RandomEffect = c(intercepts, slopes),
          Type = rep(c("Intercept", "Slope"), each = length(intercepts)),
          GroupFactor = rep(random_effect, length(intercepts) * 2)
        )
      } else {
        # If no slopes, just use intercepts
        ranef_df <- data.frame(
          Group = rownames(random_effect_data),
          RandomEffect = random_effect_data[, "(Intercept)"],
          Type = "Intercept",
          GroupFactor = random_effect
        )
      }
      
      # Merge random effects with the original data
      merged_data <- merge(data, ranef_df, by.x = random_effect, by.y = "Group", all.x = TRUE)
      
      # Convert the random effect grouping variable into a factor for color aesthetic
      merged_data[[random_effect]] <- as.factor(merged_data[[random_effect]])
      
      # If plot_response is TRUE, use the response variable (Y) for the y-axis
      if (plot_response) {
        p <- ggplot(merged_data, aes_string(x = fixed_effects, y = "Y", color = random_effect)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "black") +  # Add a linear regression line
          geom_line(aes(y = Predicted), linetype = "dashed") +  # Add model-predicted regression line
          labs(
            title = paste("Random Effects vs Response Variable (Y)"),
            x = paste(fixed_effects, "(Fixed Effect)"),
            y = "Response Variable (Y)",
            color = paste("Levels of", random_effect)
          ) +
          theme_minimal() +
          theme(legend.position = "bottom")  # Move legend to the bottom
        
        if (facet) {
          facet_labels <- paste("Group", as.numeric(unique(merged_data[[random_effect]])))
          p <- p + facet_wrap(~ merged_data[[random_effect]], labeller = labeller(.default = setNames(facet_labels, levels(merged_data[[random_effect]]))))
        }
        
        print(p)
        
      } else {
        # Plot against random effect
        p <- ggplot(merged_data, aes_string(x = fixed_effects, y = "RandomEffect", color = random_effect)) +
          geom_point() + 
          geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "black") +
          geom_line(aes(y = Predicted), linetype = "dashed") +  # Add model-predicted regression line
          labs(
            title = paste("Random Effects vs", fixed_effects),
            x = paste(fixed_effects, "(Fixed Effect)"),
            y = "Random Effect",
            color = paste("Levels of", random_effect)
          ) +
          theme_minimal() +
          theme(legend.position = "bottom")  # Move legend to the bottom
        
        if (facet) {
          facet_labels <- paste("Group", as.numeric(unique(merged_data[[random_effect]])))
          p <- p + facet_wrap(~ merged_data[[random_effect]], labeller = labeller(.default = setNames(facet_labels, levels(merged_data[[random_effect]]))))
        }
        
        print(p)
      }
    } else {
      stop("Invalid random_effect or random_effect is NULL.")
    }
  }
}

# Example usage with a GLMM model:
# Fit the GLMM model with random intercepts and slopes
PoisMod <- glmer(Y ~ X1 + X2 + X3 + (1 | Z1) + (1 | Z2) + (X1 | Z3), 
                 data = sim.data3, 
                 family = poisson(link = "log"))

# Call the plotting function for the fixed effect "X2" and random intercept for "Z2"
plot_random_effects_vs_fixed(PoisMod, sim.data3, 
                             fixed_effects = "X2", random_effect = "Z2", 
                             include_slopes = FALSE)


# Call the plotting function for the fixed effect "X1" and random intercept for "Z1", but plot against the response variable (Y)
plot_random_effects_vs_fixed(PoisMod, sim.data3, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE)


# Call the plotting function for the fixed effect "X2" and random intercept for "Z1", with faceting
plot_random_effects_vs_fixed(PoisMod, sim.data3, 
                             fixed_effects = "X2", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, 
                             facet = TRUE)

# Example fixed effects for which you want to plot the summary statistics
combined_plot <- plot_random_effects_vs_fixed(PoisMod, 
                                              sim.data3, 
                                              fixed_effects = c("X1", "X2"), # Define 1 to the max number of fixed effects
                                              random_effect = "Z1", 
                                              plot_summary_stats = TRUE, 
                                              facet = TRUE)
print(combined_plot)


# No option for 0 random effects - so GLM data