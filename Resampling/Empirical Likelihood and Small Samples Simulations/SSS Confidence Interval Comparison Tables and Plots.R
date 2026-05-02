### Small Sample Simulations - Confidence Interval Comparison Tables and Plots
# Load simulation objects
for (i in 1:3) {
  # So not to overload current environment?
  e <- new.env()
  load(paste0("./Resampling Code//SSS", i, ".RData"), envir = e)
  
  # Rename because the name of the saved objects are the same...
  if (i == 1){
    array_one <- e$int_array
  } else if (i == 2) {
    array_two <- e$int_array
  } else if (i == 3) {
    array_three <- e$int_array
  } else {
    stop()
  }
  # Cleanup
  rm(e)
  gc()
}

# Combine arrays
library(abind)
int_array <- abind(array_one, array_two, array_three)

# Remove stand alone arrays from environment for space reasons
objs <- c("array_one", "array_two", "array_three")
rm(list = objs)

### Related Simulation Information
alpha <- 0.05

### Across distributions and sample sizes, comparing methods:
# At a given target nominal alpha, matching observed coverages to compare lengths
match_obs_cov_comp_lengths <- function(results_array, alph) {
  # Assumes array is set up as:
  # nsim x alpha grid x lower/upper/captured/length/alpha x distributions x methods x sample sizes
  array_dim <- dim(results_array)
  nsim_length <- array_dim[1]
  alpha_grid_length <- array_dim[2]
  values_length <- array_dim[3]
  dists_length <- array_dim[4]
  methods_length <- array_dim[5]
  ss_length <- array_dim[6]
  
  # Function not general, so shouldn't be used out of this context
  if (is.null(dimnames(results_array))) {
    stop("Function relies on specific naming conventions")
  }
  
  # Average over simulations
  mean_results_array <- colMeans(results_array)
  
  # Fabricate output object
  output_array <- array(NA, dim = c(methods_length, methods_length, dists_length, ss_length))
  for (s in 1:ss_length) {
    for (d in 1:dists_length) {
      for (mr in 1:methods_length) {
        # Find the length and nominal alpha where our observed coverage is closest to our target
        ref <- mean_results_array[, , d, mr, s]
        closest_ref <- ref[which.min(abs(ref[, "Cap?"] - (1 - alph))), ]
        
        for (mc in 1:methods_length) {
          # Skip if the reference and comparison methods are the same
          if (mr == mc) {
            # Store length of reference method
            output_array[mc, mr, d, s] <- closest_ref["Length"]
            next
          }
          
          # Find length of the other methods at the closest observed coverage to the reference
          comp <- mean_results_array[, , d, mc, s]
          output_array[mc, mr, d, s] <- comp[which.min(abs(comp[, "Cap?"] - closest_ref["Cap?"])), ]["Length"]
        }
      }
    }
  }
  
  # Naming
  dimension_names <- dimnames(results_array)
  method_names <- dimension_names[[5]]
  dist_names <- dimension_names[[4]]
  sample_sizes <- dimension_names[[6]]
  dimnames(output_array) <- list(method_names, method_names, dist_names, sample_sizes)
  
  # Print latex code for table(s)
  
  
  return(output_array)
}

# This doesn't work well for very small sample sizes (3, 4, 5) since the highest coverage they attain is around 0.80
# should expand alpha grid and rerun if this is a problem
match_obs_cov_comp_lengths(int_array, alph = alpha)

# Average length of method when it achieves 1 - alph coverage
length_at_nominal <- function(results_array, alph) {
  # Assumes array is set up as:
  # nsim x alpha grid x lower/upper/captured/length/alpha x distributions x methods x sample sizes
  array_dim <- dim(results_array)
  nsim_length <- array_dim[1]
  alpha_grid_length <- array_dim[2]
  values_length <- array_dim[3]
  dists_length <- array_dim[4]
  methods_length <- array_dim[5]
  ss_length <- array_dim[6]
  
  # Function not general, so shouldn't be used out of this context
  if (is.null(dimnames(results_array))) {
    stop("Function relies on specific naming conventions")
  }
  
  # Average over simulations
  mean_results_array <- colMeans(results_array)
  
  # Fabricate output object
  output_array <- array(NA, dim = c(methods_length, dists_length, ss_length))
  for (s in 1:ss_length) {
    for (d in 1:dists_length) {
      for (m in 1:methods_length) {
        # Find the length and nominal alpha where our observed coverage is closest to our target
        ref <- mean_results_array[, , d, m, s]
        closest_ref <- ref[which.min(abs(ref[, "Cap?"] - (1 - alph))), ]
        
        # Store length of method
        output_array[m, d, s] <- closest_ref["Length"]
      }
    }
  }
  
  # Naming
  dimension_names <- dimnames(results_array)
  method_names <- dimension_names[[5]]
  dist_names <- dimension_names[[4]]
  sample_sizes <- dimension_names[[6]]
  dimnames(output_array) <- list(method_names, dist_names, sample_sizes)
  
  # Print latex code for table(s)
  
  
  return(output_array)
}
length_at_nominal(int_array, alph = alpha)

# Matching lengths and comparing coverage (on a dataset level /  averaged and nominal / calibrated alpha)
match_length_comp_cov <- function(results_array, alph, dataset_level = TRUE, nominal_alpha = TRUE) {
  # Assumes array is set up as:
  # nsim x alpha grid x lower/upper/captured/length/alpha x distributions x methods x sample sizes
  dimension_names <- dimnames(results_array)
  array_dim <- dim(results_array)
  nsim_length <- array_dim[1]
  alpha_grid_length <- array_dim[2]
  values_length <- array_dim[3]
  dists_length <- array_dim[4]
  methods_length <- array_dim[5]
  ss_length <- array_dim[6]
  
  # Function not general, so shouldn't be used out of this context
  if (is.null(dimnames(results_array))) {
    stop("Function relies on specific naming conventions")
  }
  
  # Dataset level requires that nominal alpha is observed
  if (!(paste("Alpha", alph) %in% dimension_names[[2]])) {
    stop("Nominal alpha specified not in alpha grid")
    # I should be able to just find the next closest alpha but I'm having a stroke trying to code it
    # to not be dependent on the naming of the array in this respect
  }
  
  # Averaged over simulation array object
  mean_int_array <- colMeans(results_array)
  
  if (dataset_level == TRUE && nominal_alpha == TRUE) {
    ### Dataset level comparison at a nominal alpha level
    # Fabricate output object
    output_array <- array(NA, dim = c(methods_length, methods_length, dists_length, ss_length))
    for (s in 1:ss_length) {
      for (d in 1:dists_length) {
        for (mr in 1:methods_length) {
          # Our reference data the nominal alpha 
          ref_data <- results_array[, paste("Alpha", alph), , d, mr, s]
          
          for (mc in 1:methods_length) {
            # Our comparison data the nominal alpha 
            comp_data <- results_array[, paste("Alpha", alph), , d, mc, s]
            
            # Temporary vector for captured indicators for a given comparison
            cap_vec <- vector("numeric", length = nsim_length)
            for (i in 1:nsim_length) {
              # Store capture indicator for comparison method at the closest length to reference
              cap_vec[i] <- comp_data[which.min(abs(comp_data[, "Length"] - ref_data[i, "Length"])), ]["Cap?"]
            }
            # Store coverage (average of indicator) for comparison method from matching the ith dataset length
            output_array[mc, mr, d, s] <- mean(cap_vec)
          }
        }
      }
    }
    
  } else if (dataset_level == TRUE && nominal_alpha == FALSE) {
    ### Dataset level comparison at a calibrated alpha level
    # Fabricate output object
    output_array <- array(NA, dim = c(methods_length, methods_length, dists_length, ss_length))
    for (s in 1:ss_length) {
      for (d in 1:dists_length) {
        for (mr in 1:methods_length) {
          # What is the nominal alpha in which we achieve our target coverage?
          temp_obj <- mean_int_array[, , d, mr, s]
          target_alpha <- temp_obj[(which.min(abs(temp_obj[, "Cap?"] - (1 - alph)))), ]["Alpha"]
          
          # Our reference data the target alpha
          ref_data <- results_array[, paste("Alpha", round(target_alpha, 5)), , d, mr, s]
          
          for (mc in 1:methods_length) {
            # Our comparison data the target alpha 
            comp_data <- results_array[, paste("Alpha", round(target_alpha, 5)), , d, mc, s]
            
            # Temporary vector for captured indicators for a given comparison
            cap_vec <- vector("numeric", length = nsim_length)
            for (i in 1:nsim_length) {
              # Store capture indicator for comparison method at the closest length to reference
              cap_vec[i] <- comp_data[which.min(abs(comp_data[, "Length"] - ref_data[i, "Length"])), ]["Cap?"]
            }
            # Store coverage (average of indicator) for comparison method from matching the ith dataset length
            output_array[mc, mr, d, s] <- mean(cap_vec)
          }
        }
      }
    }
    
  } else if (dataset_level == FALSE && nominal_alpha == TRUE) {
    ### Averaged over datasets level comparison at a nominal alpha level
    # Fabricate output object
    output_array <- array(NA, dim = c(methods_length, methods_length, dists_length, ss_length))
    for (s in 1:ss_length) {
      for (d in 1:dists_length) {
        for (mr in 1:methods_length) {
          # Reference data for the nominal alpha
          temp_obj <- mean_int_array[, , d, mr, s]
          ref_data <- temp_obj[(which.min(abs(temp_obj[, "Alpha"] - alph))), ]
          
          for (mc in 1:methods_length) {
            # Our comparison data the target alpha
            Temp_obj <- mean_int_array[, , d, mc, s]
            # Equating average length between methods what is the observed coverage rate?
            comp_data <- Temp_obj[(which.min(abs(Temp_obj[, "Length"] - ref_data["Length"]))), ]
            # Store coverage
            output_array[mc, mr, d, s] <- comp_data["Cap?"]
          }
        }
      }
    }
    
  } else if (dataset_level == FALSE && nominal_alpha == FALSE) {
    ### Averaged over datasets level comparison at a calibrated alpha level
    # Fabricate output object
    output_array <- array(NA, dim = c(methods_length, methods_length, dists_length, ss_length))
    for (s in 1:ss_length) {
      for (d in 1:dists_length) {
        for (mr in 1:methods_length) {
          # Reference data for the nominal alpha in which observed coverage is closest to target
          temp_obj <- mean_int_array[, , d, mr, s]
          ref_data <- temp_obj[(which.min(abs(temp_obj[, "Cap?"] - (1 - alph)))), ]
          
          for (mc in 1:methods_length) {
            # Our comparison data the target alpha
            Temp_obj <- mean_int_array[, , d, mc, s]
            # Equating average length between methods what is the observed coverage rate?
            comp_data <- Temp_obj[(which.min(abs(Temp_obj[, "Length"] - ref_data["Length"]))), ]
            # Store coverage
            output_array[mc, mr, d, s] <- comp_data["Cap?"]
          }
        }
      }
    }
    
  } else {
    stop("Some has gone awry!")
  }
  
  # Naming
  method_names <- dimension_names[[5]]
  dist_names <- dimension_names[[4]]
  sample_sizes <- dimension_names[[6]]
  dimnames(output_array) <- list(method_names, method_names, dist_names, sample_sizes)
  
  # Print latex code for table(s)
  
  
  return(output_array)
}
match_length_comp_cov(int_array, alph = alpha, dataset_level = FALSE, nominal_alpha = FALSE)

#### Plotting
library(ggplot2)
library(reshape2)

# Average over simulations
mean_int_array <- colMeans(int_array)

# Setting up the plotting dataframe
df_long <- melt(mean_int_array,
                varnames = c("Alpha", "Metric", "Distribution", "Method", "SampleSize"),
                value.name = "Value")
df_long[["Alpha"]] <- as.numeric(sub("Alpha", "", df_long[["Alpha"]]))
df_long[["SampleSize"]] <- as.numeric(sub("SS", "", df_long[["SampleSize"]]))

# What is the observed coverage or length for each method under each distribution at a given alpha
smol_sim_plot <- function(df, yvar, alph) {
  plot_data <- df[abs(df[["Alpha"]] - alph) < 1e-10 & df_long[["Metric"]] == yvar, ]
  
  # Base plot
  p <- ggplot(plot_data, aes(x = SampleSize, y = Value, color = Method, group = Method)) +
    geom_line(size = 1) +
    geom_point(size = 2, alpha = 0.50) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~ Distribution, nrow = 2, ncol = 4) +
    labs(
      x = "Sample Size",
      y = yvar,
      color = "Method"
    ) +
    coord_cartesian(ylim = c(0, 5)) + # Zoom without excluding points
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.key.size = unit(0.80, "cm"),
      legend.spacing.x = unit(0.50, "cm"),
      legend.spacing.y = unit(0.50, "cm"),
      legend.background = element_rect(fill = "white", color = "black")
    )
  
  # Conditionally add nominal coverage line
  if (yvar == "Cap?") {
    p <- p + geom_hline(yintercept = 1 - alph, linetype = "dashed", color = "black") +
      coord_cartesian(ylim = c(0, 1))
  }
  
  return(p)
}

# Plot Coverage
smol_sim_plot(df = df_long, yvar = "Cap?", alph = alpha)

# Plot Length
smol_sim_plot(df = df_long, yvar = "Length", alph = alpha)