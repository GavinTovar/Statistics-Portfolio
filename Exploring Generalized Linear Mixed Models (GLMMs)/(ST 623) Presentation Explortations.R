# More Packages
library(knitr)
library(kableExtra)

# Some Functions specific to this example
Extract.Coeffs <- function(model) {
  # Get the fixed effects and their standard errors
  fixed_effects <- fixef(model)
  fixed_se <- sqrt(diag(vcov(model)))
  
  # Combine fixed effects and standard errors into a matrix
  fixed_matrix <- cbind(fixed_effects, fixed_se)
  rownames(fixed_matrix) <- names(fixed_effects)
  colnames(fixed_matrix) <- c("Estimate", "Std. Error")
  
  # Get the random effects and their standard errors
  random_effects <- as.data.frame(VarCorr(model))
  
  # Create a matrix for random effects
  random_matrix <- NULL
  if (nrow(random_effects) > 0) {
    random_matrix <- matrix(ncol = 2, nrow = nrow(random_effects))
    random_matrix[, 1] <- random_effects$vcov
    random_matrix[, 2] <- ifelse(random_effects$vcov > 0, sqrt(random_effects$vcov), as.numeric(NA))
    rownames(random_matrix) <- paste("Random:", random_effects$grp)
    colnames(random_matrix) <- c("Variance", "Std. Dev.")
  }
  
  # Return a list containing both fixed and random coefficient matrices
  return(list(fixed = fixed_matrix, random = random_matrix))
}
OD.Fun <- function(model, data) {
  # Calculate the residual deviance
  residual_deviance <- deviance(model)
  
  # Get the number of observations
  n <- nrow(data)
  
  # Check if the model is GLM or GLMM
  if (inherits(model, "glmerMod")) {
    # GLMM model (with random effects)
    num_parameters <- length(fixef(model)) + length(unlist(VarCorr(model)))  # Fixed effects + random effects
  } else if (inherits(model, "glm")) {
    # GLM model (no random effects)
    num_parameters <- length(coef(model))  # Only fixed effects, use coef() for GLM
  } else {
    stop("Model type not supported. Please provide a GLM or GLMM model.")
  }
  
  # Calculate the degrees of freedom for residuals
  df_residual <- n - num_parameters
  
  # Calculate the overdispersion parameter
  overdispersion <- residual_deviance / df_residual
  
  return(overdispersion)
}
Table.Making.LMM <- function(model){
  # Extracting the summary of the model
  model_summary <- summary(model)
  
  # Extract fixed effects estimates, standard errors, and p-values
  fixed_effects <- data.frame(
    Estimate = model_summary$coefficients[, 1],
    Std_Error = model_summary$coefficients[, 2],
    t_value = model_summary$coefficients[, 3]
  )
  
  # Extract random effects (this is a simplified version)
  random_effects <- data.frame(
    Grouping = c("Z1 (Intercept)", "Z1 (X1)", "Z1 (X2)", "Z2 (Intercept)", "Z2 (X1)", "Z2 (X2)"),
    Std_Dev = c(attr(model_summary$varcor$Z1, "stddev"), attr(model_summary$varcor$Z2, "stddev"))
  )
  
  # Create a LaTeX table for fixed effects
  fixed_effects_table <- kable(fixed_effects, format = "latex", booktabs = TRUE, caption = "Fixed Effects Model Summary") %>%
    kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, background = "lightgray")  # Add color to the header row
  
  # Create a LaTeX table for random effects
  random_effects_table <- kable(random_effects, format = "latex", booktabs = TRUE, caption = "Random Effects Model Summary") %>%
    kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, background = "lightgray")  # Add color to the header row
  
  
  # Print the tables (for RMarkdown or Quarto use, the result will be displayed automatically)
  return(list(fixed_effects_table, random_effects_table))
}
Table.Making.GLM <- function(model){
  # Extracting the summary of the model
  model_summary <- summary(model)
  
  # Extract fixed effects estimates, standard errors, and p-values
  fixed_effects <- data.frame(
    Estimate = model_summary$coefficients[, 1],
    Std_Error = model_summary$coefficients[, 2],
    t_value = model_summary$coefficients[, 3],
    p_value = model_summary$coefficients[, 4]
  )
  
  # Create a LaTeX table for fixed effects
  fixed_effects_table <- kable(fixed_effects, format = "latex", booktabs = TRUE, caption = "Fixed Effects Model Summary") %>%
    kable_styling(latex_options = "striped", full_width = F) %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(2, background = "lightgray")
  
  # Print the tables (for RMarkdown or Quarto use, the result will be displayed automatically)
  return(fixed_effects_table)
  
}
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

### Exploration
# Standard GLMM data (no overdispersion)
set.seed(10262025)
sim.data0 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(5, 4),
  N = 240,
  Sigma.int = c(0.05, 0.05),
  Sigma.slope = c(0.15, 0.15),
  Cov.mat = NULL,
  Beta = c(0.15, 0.25, 0.50),
  mu.vec = c(1, 1),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)
sim.data0

# Comparing Fits
PoisMod0 <- glmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
                  data = sim.data0, 
                  family = poisson(link = "log"))


# With GLMM data
# What happens when we fit a LMM instead (Missing a link function)
PoisLMM <- lmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
                data = sim.data0)

# Model summary
summary(PoisLMM)

# Estimates are completely off but like duh

# With GLMM data
# What happens when we fit a GLM instead (Missing random effects)
PoisGLM <- glm(Y ~ X1 + X2, data = sim.data0, family = poisson(link = "log"))

# Model summary
summary(PoisGLM)
OD.Fun(model = PoisGLM, data = sim.data0)

# it's actually kinda fine...? Not crazy overdispersion??
# Maybe because the random effects are not high in magnitude...?

# Increasing magnitude of random effects to try again
set.seed(10262025)
sim.data00 <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(5, 4),
  N = 240,
  Sigma.int = c(1, 1),
  Sigma.slope = c(0.75, 0.75),
  Cov.mat = NULL,
  Beta = c(0.15, 0.25, 0.50),
  mu.vec = c(1, 1),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = NULL
)
sim.data00

PoisLMM2 <- lmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
                 data = sim.data00)
PoisGLM2 <- glm(Y ~ X1 + X2, data = sim.data00, family = poisson(link = "log"))

summary(PoisLMM2)
summary(PoisGLM2)
OD.Fun(model = PoisGLM2, data = sim.data00)

# Add some Gaussian noise to see if idk
set.seed(10262025)
sim.data.noise <- Mixed.Model.DataGen(
  FixEff = 3,
  RanEff = 2,
  Levels = c(5, 4),
  N = 240,
  Sigma.int = c(1, 1),
  Sigma.slope = c(0.75, 0.75),
  Cov.mat = NULL,
  Beta = c(0.15, 0.25, 0.50),
  mu.vec = c(1, 1),
  Response.distribution = "Poisson",
  Error.params = NULL,
  Link.function = "log",
  overdispersion = FALSE,
  Noise = rnorm(240, mean = 0, sd = 5)
)

PoisLMMNoise <- lmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
                     data = sim.data.noise)
PoisGLMNoise <- glm(Y ~ X1 + X2, data = sim.data.noise, family = poisson(link = "log"))

summary(PoisLMMNoise)
summary(PoisGLMNoise)
OD.Fun(model = PoisGLMNoise, data = sim.data.noise)

# Yes, now that the random effect coefficients are larger - the overdispersion is higher
# Create plots to follow these explorations

# Tables for prezi
Table.Making.LMM(PoisLMM)
Table.Making.GLM(PoisGLM)
Table.Making.LMM(PoisLMM2)
Table.Making.GLM(PoisGLM2)
Table.Making.GLM(PoisMod0)

# Plotting

# Call the plotting function for the fixed effect "X2" and random intercept for "Z1", with faceting
plot_random_effects_vs_fixed(PoisMod0, sim.data0, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, 
                             facet = TRUE) + theme(
                               legend.position = c(0.90, 0.10),
                               legend.justification = c("right", "bottom"),
                               legend.title = element_text(size = 15),
                               legend.text = element_text(size = 14),
                               legend.key.size = unit(1.5, "lines")
                             )

plot_random_effects_vs_fixed(PoisLMM, sim.data0, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, 
                             facet = TRUE)

plot_random_effects_vs_fixed(glmer(Y ~ X1 + X2 + (X1 + X2 | Z1) + (X1 + X2 | Z2), 
                                   data = sim.data00, 
                                   family = poisson(link = "log")), sim.data00, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, 
                             facet = TRUE)

plot_random_effects_vs_fixed(PoisLMM2, sim.data00, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, 
                             facet = TRUE)


# Example fixed effects for which you want to plot the summary statistics
combined_plot <- plot_random_effects_vs_fixed(PoisMod0, sim.data0,
                                              fixed_effects = c("X1", "X2"), # Define 1 to the max number of fixed effects
                                              random_effect = "Z1", 
                                              plot_summary_stats = TRUE, 
                                              facet = TRUE)
print(combined_plot)




Table.Making.LMM <- function(model){
  # Extracting the summary of the model
  model_summary <- summary(PoisLMM)
  
  # Extract fixed effects estimates, standard errors, and p-values
  fixed_effects <- data.frame(
    Estimate = model_summary$coefficients[, 1],
    Std_Error = model_summary$coefficients[, 2],
    t_value = model_summary$coefficients[, 3]
  )
  
  # Extract random effects (this is a simplified version)
  random_effects <- data.frame(
    Grouping = c("Z1 (Intercept)", "Z1 (X1)", "Z1 (X2)", "Z2 (Intercept)", "Z2 (X1)", "Z2 (X2)"),
    Std_Dev = c(attr(model_summary$varcor$Z1, "stddev"), attr(model_summary$varcor$Z2, "stddev"))
  )
  
  # Create a LaTeX table for fixed effects
  fixed_effects_table <- kable(fixed_effects, format = "latex", booktabs = TRUE, caption = "Fixed Effects Model Summary") %>%
    kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, background = "lightgray")  # Add color to the header row
  
  # Create a LaTeX table for random effects
  random_effects_table <- kable(random_effects, format = "latex", booktabs = TRUE, caption = "Random Effects Model Summary") %>%
    kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, background = "lightgray")  # Add color to the header row
  
  
  # Print the tables (for RMarkdown or Quarto use, the result will be displayed automatically)
  return(list(fixed_effects_table, random_effects_table))
}
