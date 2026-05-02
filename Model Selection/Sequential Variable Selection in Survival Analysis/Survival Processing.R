# Load in data
load("./Survival Project/Survival Simulation Results/SurvSimResultsRData") # Not on GitHub

library(tidyverse)

# Levels
Correlation <- c("Independent", "ModerateCorrelation")
Censoring <- c("LowCens", "HighCens")
TrueDist <- c("weibull", "exponential", "lognormal", "cox")
FitDist <- c("AFT1_Weibull", "AFT2_Exponential", "AFT3_Lognormal")

### Helper Functions
# Function to flatten the nested list into a data frame
flatten_nested_list <- function(nested_list) {
  df_list <- list()
  
  for (corr_level in names(nested_list)) {
    for (cens_level in names(nested_list[[corr_level]])) {
      for (dist in names(nested_list[[corr_level]][[cens_level]])) {
        for (model in names(nested_list[[corr_level]][[cens_level]][[dist]])) {
          coefs <- nested_list[[corr_level]][[cens_level]][[dist]][[model]]
          
          df_list[[length(df_list) + 1]] <- data.frame(
            Correlation = corr_level,
            Censoring = cens_level,
            Distribution = dist,
            Model = model,
            t(as.data.frame(coefs))
          )
        }
      }
    }
  }
  
  # Combine all rows into a single data frame
  final_df <- bind_rows(df_list)
  return(final_df)
}

# Function to reorder columns by prefix
reorder_columns <- function(df, prefix_length = 5) {
  # Identify columns with the same prefix
  prefix_groups <- split(names(df), substr(names(df), 1, prefix_length))
  
  # Flatten the list of columns and reorder them
  new_column_order <- unlist(prefix_groups)
  
  # Reorder the columns
  df <- df[, new_column_order]
  
  return(df)
}

# Function to process and format the data frame from multiple nested lists
process_data_frame <- function(...) {
  # Combine multiple nested lists into a list
  nested_lists <- list(...)
  
  # Initialize an empty list to hold all the processed data frames
  df_list <- list()
  
  # Loop over each nested list and process it
  for (nested_list in nested_lists) {
    # Flatten the nested list
    flat_df <- flatten_nested_list(nested_list)
    
    # Renaming the 'Model' column values
    processed_df <- flat_df %>%
      mutate(Model = recode(Model,
                            "AFT1_Weibull" = "Weibull",
                            "AFT2_Exponential" = "Exponential",
                            "AFT3_Lognormal" = "Lognormal")) %>%
      arrange(Distribution, Censoring, Correlation) %>%
      group_by(Distribution, Censoring, Correlation) %>%
      mutate(Distribution = factor(Distribution, levels = c("weibull", "exponential", "lognormal", "cox"))) %>%
      select(Distribution, Censoring, Correlation, Model, everything()) %>%
      arrange(Distribution)
    
    # Append the processed data frame to the df_list
    df_list[[length(df_list) + 1]] <- processed_df
  }
  
  # If there are multiple data frames, combine them by columns
  if (length(df_list) > 1) {
    # Combine all processed data frames by column (cbind)
    combined_df <- do.call(cbind, df_list)
    
    # Select columns that start with "Intercept"
    intercept_columns <- select(combined_df, starts_with("Intercept"))
    
    # Select columns that start with "Beta"
    beta_columns <- select(combined_df, starts_with("Beta"))
    
    # Combine the info columns first, 
    # then add the "Intercept" columns after the 4th column, 
    # then the "Beta" columns, reordered by Beta type.
    final_df <- cbind(select(combined_df, 1:4), intercept_columns, reorder_columns(beta_columns, prefix_length = 5))
    
  } else {
    # If there's only one data frame, just return the processed one
    final_df <- df_list[[1]]
  }
  
  return(as.data.frame(final_df))
}

### Tables 1-5
# Initialize the lists before populating them
FullModels <- Scen11.Models <- Scen12.Models <- Scen21.Models <- Scen22.Models <- list()

# Nested loops to iterate over all combinations
for (correlation in Correlation) {
  # Initialize nested list if it doesn't exist
  if (is.null(FullModels[[correlation]])) {
    FullModels[[correlation]] <- list()
  }
  if (is.null(Scen11.Models[[correlation]])) {
    Scen11.Models[[correlation]] <- list()
  }
  if (is.null(Scen12.Models[[correlation]])) {
    Scen12.Models[[correlation]] <- list()
  }
  if (is.null(Scen21.Models[[correlation]])) {
    Scen21.Models[[correlation]] <- list()
  }
  if (is.null(Scen22.Models[[correlation]])) {
    Scen22.Models[[correlation]] <- list()
  }
  
  for (censoring in Censoring) {
    # Initialize nested list if it doesn't exist
    if (is.null(FullModels[[correlation]][[censoring]])) {
      FullModels[[correlation]][[censoring]] <- list()
    }
    if (is.null(Scen11.Models[[correlation]][[censoring]])) {
      Scen11.Models[[correlation]][[censoring]] <- list()
    }
    if (is.null(Scen12.Models[[correlation]][[censoring]])) {
      Scen12.Models[[correlation]][[censoring]] <- list()
    }
    if (is.null(Scen21.Models[[correlation]][[censoring]])) {
      Scen21.Models[[correlation]][[censoring]] <- list()
    }
    if (is.null(Scen22.Models[[correlation]][[censoring]])) {
      Scen22.Models[[correlation]][[censoring]] <- list()
    }
    
    for (truedist in TrueDist) {
      # Initialize nested list if it doesn't exist
      if (is.null(FullModels[[correlation]][[censoring]][[truedist]])) {
        FullModels[[correlation]][[censoring]][[truedist]] <- list()
      }
      if (is.null(Scen11.Models[[correlation]][[censoring]][[truedist]])) {
        Scen11.Models[[correlation]][[censoring]][[truedist]] <- list()
      }
      if (is.null(Scen12.Models[[correlation]][[censoring]][[truedist]])) {
        Scen12.Models[[correlation]][[censoring]][[truedist]] <- list()
      }
      if (is.null(Scen21.Models[[correlation]][[censoring]][[truedist]])) {
        Scen21.Models[[correlation]][[censoring]][[truedist]] <- list()
      }
      if (is.null(Scen22.Models[[correlation]][[censoring]][[truedist]])) {
        Scen22.Models[[correlation]][[censoring]][[truedist]] <- list()
      }
      
      for (fitdist in FitDist) {
        # Extract the coefficients and store them in the lists
        FullModels[[correlation]][[censoring]][[truedist]][[fitdist]] <- SimulationResults[[correlation]][[censoring]][[truedist]][[fitdist]][["Full_Coef"]] # S00
        Scen11.Models[[correlation]][[censoring]][[truedist]][[fitdist]] <- SimulationResults[[correlation]][[censoring]][[truedist]][[fitdist]][["MS_Coef"]] # S11
        Scen12.Models[[correlation]][[censoring]][[truedist]][[fitdist]] <- SimulationResults[[correlation]][[censoring]][[truedist]][[fitdist]][["MS_Sel_Coef"]] # S12 
        Scen21.Models[[correlation]][[censoring]][[truedist]][[fitdist]] <- SimulationResults[[correlation]][[censoring]][[truedist]][[fitdist]][["Sel_Coef"]] # S21
        Scen22.Models[[correlation]][[censoring]][[truedist]][[fitdist]] <- SimulationResults[[correlation]][[censoring]][[truedist]][[fitdist]][["Sel_MS_Coef"]] # S22
      } 
    }
  }
}

# After filling the lists, you can now compute column means for each combination
# Table 1 - Full models - Scenario 00
column_means_full <- list()
for (correlation in Correlation) {
  for (censoring in Censoring) {
    for (truedist in TrueDist) {
      for (fitdist in FitDist) {
        # Ensure the model exists
        if (!is.null(FullModels[[correlation]][[censoring]][[truedist]][[fitdist]])) {
          coef_matrix <- FullModels[[correlation]][[censoring]][[truedist]][[fitdist]]
          # Calculate column means of the coefficient matrix (na.rm = TRUE to ignore NA values)
          column_means_full[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(coef_matrix, na.rm = TRUE)
        }
      }
    }
  }
}

# Column means 'FullModels'
column_means_full

# Table 2 - Scenario 11
column_means_scen11a <- column_means_scen11b <- column_means_scen11c <- list()
for (correlation in Correlation) {
  for (censoring in Censoring) {
    for (truedist in TrueDist) {
      for (fitdist in FitDist) {
        # Ensure the model exists
        if (!is.null(Scen11.Models[[correlation]][[censoring]][[truedist]][[fitdist]])) {
          coef_matrix <- Scen11.Models[[correlation]][[censoring]][[truedist]][[fitdist]]
          # Calculate column means of the coefficient matrix (na.rm = TRUE to ignore NA values)
          column_means_scen11a[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(coef_matrix, na.rm = TRUE)
          column_means_scen11b[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(replace(coef_matrix, is.na(coef_matrix), 0))
          column_means_scen11c[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(!is.na(coef_matrix))
        }
      }
    }
  }
}

# Table 3 - Scenario 12
column_means_scen12a <- column_means_scen12b <- column_means_scen12c <- list()
for (correlation in Correlation) {
  for (censoring in Censoring) {
    for (truedist in TrueDist) {
      for (fitdist in FitDist) {
        # Ensure the model exists
        if (!is.null(Scen12.Models[[correlation]][[censoring]][[truedist]][[fitdist]])) {
          coef_matrix <- Scen12.Models[[correlation]][[censoring]][[truedist]][[fitdist]]
          # Calculate column means of the coefficient matrix (na.rm = TRUE to ignore NA values)
          column_means_scen12a[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(coef_matrix, na.rm = TRUE)
          column_means_scen12b[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(replace(coef_matrix, is.na(coef_matrix), 0))
          column_means_scen12c[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(!is.na(coef_matrix))
        }
      }
    }
  }
}

# Table 4 - Scenario 21
column_means_scen21a <- column_means_scen21b <- column_means_scen21c <- list()
for (correlation in Correlation) {
  for (censoring in Censoring) {
    for (truedist in TrueDist) {
      for (fitdist in FitDist) {
        # Ensure the model exists
        if (!is.null(Scen21.Models[[correlation]][[censoring]][[truedist]][[fitdist]])) {
          coef_matrix <- Scen21.Models[[correlation]][[censoring]][[truedist]][[fitdist]]
          # Calculate column means of the coefficient matrix (na.rm = TRUE to ignore NA values)
          column_means_scen21a[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(coef_matrix, na.rm = TRUE)
          column_means_scen21b[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(replace(coef_matrix, is.na(coef_matrix), 0))
          column_means_scen21c[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(!is.na(coef_matrix))
        }
      }
    }
  }
}

# Table 5 - Scenario 22
column_means_scen22a <- column_means_scen22b <- column_means_scen22c <- list()
for (correlation in Correlation) {
  for (censoring in Censoring) {
    for (truedist in TrueDist) {
      for (fitdist in FitDist) {
        # Ensure the model exists
        if (!is.null(Scen22.Models[[correlation]][[censoring]][[truedist]][[fitdist]])) {
          coef_matrix <- Scen22.Models[[correlation]][[censoring]][[truedist]][[fitdist]]
          # Calculate column means of the coefficient matrix (na.rm = TRUE to ignore NA values)
          column_means_scen22a[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(coef_matrix, na.rm = TRUE)
          column_means_scen22b[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(replace(coef_matrix, is.na(coef_matrix), 0))
          column_means_scen22c[[correlation]][[censoring]][[truedist]][[fitdist]] <- colMeans(!is.na(coef_matrix))
        }
      }
    }
  }
}


# Apply the function to different input lists
Table1 <- process_data_frame(column_means_full)
Table2 <- process_data_frame(column_means_scen11a)
Table3 <- process_data_frame(column_means_scen12a, column_means_scen12b, column_means_scen12c)
Table4 <- process_data_frame(column_means_scen21a, column_means_scen21b, column_means_scen21c)
Table5 <- process_data_frame(column_means_scen22a, column_means_scen22b, column_means_scen22c)

require(xtable)
### xtable for Table 1 or 2 -- 
df <- data.frame(
  Distribution=rep(c(""), 48),
  Censoring=rep(c(""), 48),
  Correlation=rep(c(""), 24),
  Model=rep(c("Weibull", "Exponential", "Lognormal"), 16),
  Coefficients=round(Table2[,-(1:4)], 3)
)

row.lengths.1 <- rep(12, 4)
df[[1]][1+(0:3)*12] <- paste0("\\midrule\\multirow{", row.lengths.1, "}{*}{\\textbf{", c("Weib", "Exp", "Lognorm", "Gomp"), "}}")

row.lengths.2 <- rep(6, 8)
df[[2]][1+(0:7)*6] <- paste0("\\multirow{", row.lengths.2, "}{*}{\\textbf{", rep(c("High", "Low"), 4), "}}")

row.lengths.3 <- rep(3, 16)
df[[3]][1+(0:15)*3] <- paste0("\\multirow{", row.lengths.3, "}{*}{\\textbf{", rep(c("Ind", "Cor"), 6), "}}")

strCaption <- paste0("\\textbf{Table Whatever} This table is just produced with some ",
                     "random data and does not mean anything. Just to show you how ",
                     "things work.")

print(xtable(df, digits = c(0, 0, 0, 0, rep(3, 6)), # first zero "represents" row numbers which we skip later
             align = "lllll|rrrrr",  # align and put a vertical line (first "l" again represents column of row numbers)
             caption = strCaption, label = paste0("Table", 3)),
      size = "footnotesize", #Change size; useful for bigger tables "normalsize" "footnotesize"
      include.rownames = FALSE, #Don't print rownames
      include.colnames = FALSE, #We create them ourselves
      caption.placement = "top", #"top", NULL
      hline.after=NULL, #We don't need hline; we use booktabs
      floating=TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
      sanitize.text.function = force, # Important to treat content of first column as latex function
      add.to.row = list(pos = list(-1, nrow(df)),
                        command = c(paste("\\toprule \n",  
                                          "Dist & Cens & Cor & Model & $\\hat{\\beta}_0$ & $\\hat{\\beta}_1$ & $\\hat{\\beta}_2$ & $\\hat{\\beta}_3$ & $\\hat{\\beta}_4$ \\\\ \n", 
                                          "\\midrule \n"),
                                    paste("\\bottomrule \n")  # paste is used as it is more flexible regarding adding lines
                                    )
                        )
)


### xtable for Table 3-5
df <- data.frame(
  Distribution = rep(c(""), 48),
  Censoring = rep(c(""), 48),
  Correlation = rep(c(""), 24),
  Model = rep(c("Weibull", "Exponential", "Lognormal"), 16),
  Coefficients = round(Table5[, -(1:4)], 3)
)

# Adjust percentage columns (5 sets of 3 = 15 columns)
df[, c(4 + 3 * (1:5))] <- 100 * df[, c(4 + 3 * (1:5))]

# Add Distribution multirow
row.lengths.1 <- rep(12, 4)
df[[1]][1 + (0:3) * 12] <-
  paste0("\\midrule\\multirow{", row.lengths.1, "}{*}{\\textbf{", c("Weib", "Exp", "Lognorm", "Gomp"), "}}")

# Add Censoring multirow
row.lengths.2 <- rep(6, 8)
df[[2]][1 + (0:7) * 6] <-
  paste0("\\multirow{", row.lengths.2, "}{*}{\\textbf{", rep(c("High", "Low"), 4), "}}")

# Add Correlation multirow
row.lengths.3 <- rep(3, 16)
df[[3]][1 + (0:15) * 3] <-
  paste0("\\multirow{", row.lengths.3, "}{*}{\\textbf{", rep(c("Ind", "Cor"), 6), "}}")

# Caption
strCaption <- paste0("\\textbf{Table Whatever} This table is just produced with some ",
                     "random data and does not mean anything. Just to show you how ",
                     "things work.")

# xtable print
print(xtable(df, 
             digits = c(0, 0, 0, 0, 0, rep(c(2, 2, 0), 5)),
             align = "lllll|ccc|ccc|ccc|ccc|ccc",
             caption = strCaption, 
             label = paste0("Table", 5)),
      size = "footnotesize",
      include.rownames = FALSE,
      include.colnames = FALSE,
      caption.placement = "top",
      hline.after = NULL,
      floating = TRUE,
      sanitize.text.function = force,
      add.to.row = list(
        pos = list(-1, nrow(df)),
        command = c(
          paste(
            "\\toprule \n",
            "Dist & Cens & Cor & Model & \\multicolumn{3}{c}{$\\hat{\\beta}_0 \\left(\\beta_{0} = 0.10\\right)$} & \\multicolumn{3}{c}{$\\hat{\\beta}_1 \\left(\\beta_{1} = 0.20\\right)$} & \\multicolumn{3}{c}{$\\hat{\\beta}_2 \\left(\\beta_{2} = 0.15\\right)$} 
            & \\multicolumn{3}{c}{$\\hat{\\beta}_3 \\left(\\beta_{3} = 0.075\\right)$} & \\multicolumn{3}{c}{$\\hat{\\beta}_4 \\left(\\beta_{4} = 0\\right)$} \\\\ \n",
            " & & & & NA & 0 & \\% & NA & 0 & \\% & NA & 0 & \\% & NA & 0 & \\% & NA & 0 & \\% \\\\ \n",
            "\\midrule \n"
          ),
          "\\bottomrule \n"
        )
      )
)
