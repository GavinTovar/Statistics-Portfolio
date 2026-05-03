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
        geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
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
          geom_point() +  # Plot data points
          geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +  # Add a linear regression line
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
          geom_point() +  # Plot data points
          geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +  # Add a linear regression line
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
                 data = sim.data2, 
                 family = poisson(link = "log"))


# Call the plotting function for the fixed effect "X1" and random intercept for "Z1"
plot_random_effects_vs_fixed(PoisMod, sim.data2, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE)


# Call the plotting function for the fixed effect "X1" and random intercept for "Z1", but plot against the response variable (Y)
plot_random_effects_vs_fixed(PoisMod, sim.data2, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE)


# Call the plotting function for the fixed effect "X1" and random intercept for "Z1", with faceting
plot_random_effects_vs_fixed(PoisMod, sim.data2, 
                             fixed_effects = "X1", random_effect = "Z1", 
                             include_slopes = FALSE, plot_response = TRUE, facet = TRUE)

# Example fixed effects for which you want to plot the summary statistics
fixed_effects <- c("mean_X1", "mean_X2", "mean_X3")

# Call the function to plot random effects vs fixed effects with summary statistics
fixed_effects <- c("X1", "X2")  # Define multiple fixed effects
combined_plot <- plot_random_effects_vs_fixed(PoisMod, 
                                              sim.data2, 
                                              fixed_effects = fixed_effects, 
                                              random_effect = "Z2", 
                                              plot_summary_stats = TRUE, 
                                              facet = TRUE)
print(combined_plot)
