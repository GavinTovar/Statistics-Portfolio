# Load data
# load("./Survival Project/Survival Simulation Results/SurvSimTabulated.RData")
# load("./Survival Project/Survival Simulation Results/SurvSimChosenMod.RData")

# Load some packages and functions
source("./Survival Project/R Code/Survival - Functions and Packages.R")
library(tidyverse)
library(patchwork)
library(colorspace)
library(cowplot)

# Survival Simulations - Plotting Code
SmooshedTable <- as.data.frame(bind_rows(Tabulated.Results))

# Create true coef and bias columns (assumes first p columns are the coefficients)
Beta <- c(c(0.1, 0.2), c(0.15, 0.075, 0)) 

for (i in 1:length(Beta)) {
  # True coefficient values
  SmooshedTable[[paste0("Beta_", i - 1, "_True")]] <- Beta[i]
  
  # Bias columns
  SmooshedTable[[paste0("Beta_", i - 1, "_Bias")]] <- SmooshedTable[[i]] - Beta[i]
}

# Remove Intercept columns
SmooshedTable <- select(SmooshedTable, -c(matches("^Beta_0"), matches("^Intercept")))

# Subset data into different unique scenario / simulation settings
SubPlotDataList <- SmooshedTable %>%
  group_by(Scenario, Correlation, Censoring, Gen.Distribution) %>%
  group_split()

# Marginalize over regressor variable
Smoosh.Fun <- function(df_list) {
  result_list <- vector("list", length(df_list))
  
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    
    # Split into list
    groups <- df %>%
      group_by(Distribution, Averaging) %>%
      group_split()
    
    # Average each group
    smoothed_groups <- lapply(groups, function(sub_df) {
      num_cols <- sapply(sub_df, is.numeric)
      
      # Take means of numeric columns
      averaged <- as.data.frame(lapply(sub_df[, num_cols, drop = FALSE], mean))
      
      # Take first row of non-numeric (meta) columns (removing regressor column)
      meta <- sub_df[1, !num_cols, drop = FALSE][-1]
      
      # Combine and return
      cbind(meta, averaged)
    })
    
    # Combine all regressor groups for this original data frame
    result_list[[i]] <- do.call(rbind, smoothed_groups)
  }
  
  return(result_list)
}

# Breaking into 5 lists, one for each scenario
Scenario00SubPlotList <- Smoosh.Fun(SubPlotDataList[1:16])
Scenario11SubPlotList <- Smoosh.Fun(SubPlotDataList[17:32])
Scenario12SubPlotList <- Smoosh.Fun(SubPlotDataList[33:48])
Scenario21SubPlotList <- Smoosh.Fun(SubPlotDataList[49:64])
Scenario22SubPlotList <- Smoosh.Fun(SubPlotDataList[65:80])


LinePlot.Fun <- function(PlottingData, yvar = "Bias") {
  if (yvar == "Coverage") {
    yvar <- "Captured"
  }
  
  # Create list to store individual plots
  individual_plots <- list()
  
  # Define possible levels for the distribution
  da_levels <- c("AFT1_Weibull", "AFT2_Exponential", "AFT3_Lognormal", "Cox_Regression")
  
  # Determine the range for consistent y-axis limits
  y_range <- range(unlist(lapply(PlottingData, function(df) Reshape.Fun(df)[[yvar]])), na.rm = TRUE)
  
  for (i in seq_along(PlottingData)) {
    PlotDatas <- Reshape.Fun(PlottingData[[i]])
    
    # Get present levels dynamically
    present_levels <- unique(na.omit(PlotDatas[["Distribution"]]))
    num_legends <- length(present_levels)
    
    Scen <- case_when(
      PlotDatas[["Scenario"]][1] == "S00" ~ "Full Model Fitting",
      PlotDatas[["Scenario"]][1] == "S11" ~ "Model Selection via Information Criterion",
      PlotDatas[["Scenario"]][1] == "S12" ~ "Model Selection then Variable Selection",
      PlotDatas[["Scenario"]][1] == "S21" ~ "Variable Selection via Backwards Selection",
      PlotDatas[["Scenario"]][1] == "S22" ~ "Variable Selection then Model Selection",
      TRUE ~ "Error"
    )
    
    Cens <- ifelse(PlotDatas[["Censoring"]][1] == "LowCens", "Low Cens.", "High Cens.")
    Cor <- ifelse(PlotDatas[["Correlation"]][1] == "Independent", "Indep.", "Mod. Corr.")
    Gen.Dist <- case_when(
      PlotDatas[["Gen.Distribution"]][1] == "weibull" ~ "Weibull Distribution",
      PlotDatas[["Gen.Distribution"]][1] == "exponential" ~ "Exponential Distribution",
      PlotDatas[["Gen.Distribution"]][1] == "lognormal" ~ "LogNormal Distribution",
      PlotDatas[["Gen.Distribution"]][1] == "cox" ~ "Gompertz Distribution",
      TRUE ~ "Error"
    )
    
    true_levels <- sort(unique(PlotDatas$True))
    
    # Custom x-axis label formatting
    custom_label <- function(x) {
      ifelse(x == 0, "0", formatC(x, format = "f", digits = 3, drop0trailing = TRUE))
    }
    
    # Color mapping for each yvar
    color_mapping <- if (yvar == "Bias") {
      c("AFT1_Weibull" = "#fa6d46",
        "AFT2_Exponential" = "#70e000",
        "AFT3_Lognormal" = "#00bbf9",
        "Cox_Regression" = "#9d4edd")
    } else {
      c("AFT1_Weibull" = "#c7512e",
        "AFT2_Exponential" = "#4ba600",
        "AFT3_Lognormal" = "#007ea7",
        "Cox_Regression" = "#6a1abb")
    }
    
    # Build the plot
    individual_plots[[i]] <- ggplot(PlotDatas) +
      geom_line(aes(x = True, y = !!sym(yvar), color = Distribution, linetype = Averaging), linewidth = 0.75) +
      geom_point(aes(x = True, y = !!sym(yvar), color = Distribution, shape = Averaging)) +
      {if (yvar == "Bias") geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black")} +
      {if (yvar == "Captured") geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "black")} +
      scale_x_continuous(breaks = true_levels, labels = custom_label) +
      scale_y_continuous(limits = y_range, breaks = pretty(y_range, n = 4)) +
      scale_color_manual(name = "Model", values = color_mapping[present_levels], drop = TRUE) +
      {if (num_legends < 3) guides(color = "none")} +
      labs(
        y = paste(Cens, "\n", Cor),
        x = "Data Generating \n Coefficient Value",
        subtitle = Gen.Dist
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10)
      )
    
    # Remove subtitle for plots not in top row
    if (i > 4) {
      individual_plots[[i]] <- individual_plots[[i]] + labs(subtitle = "")
    }
    
    # Remove Y-axis text for non-left plots
    if (!(i %in% c(1, 5, 9, 13))) {
      individual_plots[[i]] <- individual_plots[[i]] +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(y = "")
    }
    
    # Remove X-axis text for non-bottom plots
    if (i < 13) {
      individual_plots[[i]] <- individual_plots[[i]] +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(x = "")
    }
  }
  
  # Combine all plots into grid
  combined_plot <- wrap_plots(individual_plots, ncol = 4) + 
    plot_layout(guides = "collect") & 
    plot_annotation(
      title = paste(ifelse(yvar == "Bias", "Post-Model Selection Bias", "Confidence Interval Coverage"),
                    "for Scenario", PlottingData[[1]][["Scenario"]][1]),
      caption = Scen,
      theme = theme(
        plot.margin = margin(2, 2, 0.5, 2),
        plot.caption = element_text(size = 10, face = "italic", hjust = 0, vjust = 3),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "bottom"
      )
    )
  
  return(combined_plot)
}
# Bias Plots
pB0 <- LinePlot.Fun(Scenario00SubPlotList, yvar = "Bias")
pB1 <- LinePlot.Fun(Scenario11SubPlotList, yvar = "Bias")
pB2 <- LinePlot.Fun(Scenario12SubPlotList, yvar = "Bias")
pB3 <- LinePlot.Fun(Scenario21SubPlotList, yvar = "Bias")
pB4 <- LinePlot.Fun(Scenario22SubPlotList, yvar = "Bias")
# Coverage Plots
pC0 <- LinePlot.Fun(Scenario00SubPlotList, yvar = "Captured")
pC1 <- LinePlot.Fun(Scenario11SubPlotList, yvar = "Captured")
pC2 <- LinePlot.Fun(Scenario12SubPlotList, yvar = "Captured")
pC3 <- LinePlot.Fun(Scenario21SubPlotList, yvar = "Captured")
pC4 <- LinePlot.Fun(Scenario22SubPlotList, yvar = "Captured")

# Scenario 00
full_plot <- wrap_plots(pB0, pC0, ncol = 2) + 
  plot_annotation(
    title = "Full Model 'Bias' and Coverage by Simulation Setting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )
full_plot

# Scenario 11
scen11_plot <- wrap_plots(pB1, pC1, ncol = 2) + 
  plot_annotation(
    title = "Scenario 1 - Part 1: 'Bias' and Coverage by Simulation Setting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )
scen11_plot

# Scenario 12
scen12_plot <- wrap_plots(pB2, pC2, ncol = 2) + 
  plot_annotation(
    title = "Scenario 1 - Part 2: 'Bias' and Coverage by Simulation Setting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )
scen12_plot

# Scenario 21
scen21_plot <- wrap_plots(pB3, pC3, ncol = 2) + 
  plot_annotation(
    title = "Scenario 2 - Part 1: 'Bias' and Coverage by Simulation Setting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )
scen21_plot

# Scenario 22
scen22_plot <- wrap_plots(pB4, pC4, ncol = 2) + 
  plot_annotation(
    title = "Scenario 2 - Part 2: 'Bias' and Coverage by Simulation Setting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )
scen22_plot

### Heatmap Action
# Remove coefficient columns
SubPlotDataList2 <- lapply(SubPlotDataList, function(x) select(x, -c(matches("^Beta"), matches("^Intercept"))))

# Remove unconditional (could remove conditional rows - doesn't matter) since we only care about counts
SubPlotDataList2 <- lapply(SubPlotDataList2, function(x) x[!x[["Averaging"]] == "Unconditional", ])

# Breaking into 5 lists, one for each scenario
# Scenario00SubPlotList2 <- SubPlotDataList2[1:16] # All full models
Scenario11SubPlotList2 <- SubPlotDataList2[17:32] # All full models
Scenario12SubPlotList2 <- SubPlotDataList2[33:48]
Scenario21SubPlotList2 <- SubPlotDataList2[49:64]
Scenario22SubPlotList2 <- SubPlotDataList2[65:80]

# ---- Generate Heatmaps ----
HeatmapPlot.Fun <- function(PlottingData, ratio = 1, text_size = 3, tittle = "Default Title") {
  # Getting some overall list data
  stopifnot(length(PlottingData) == 16)
  FrequencyFrame <- bind_rows(PlottingData)
  
  # Optional: consistent model levels across all heatmaps
  all_models <- sort(unique(FrequencyFrame$Regressors))
  
  # Optional: consistent fill scale
  fill_range <- range(FrequencyFrame$Count, na.rm = TRUE)
  
  # Ensure consistent y-axis labels (using all possible 'Distribution' values)
  all_distributions <- sort(unique(FrequencyFrame$Distribution))
  
  heatmap_plots <- list()
  
  for (i in seq_along(PlottingData)) {
    PlotDatas <- PlottingData[[i]]
    
    # Descriptive labels
    Scen <- dplyr::case_when(
      PlotDatas$Scenario[1] == "S00" ~ "Scenario 0: Full Model Fitting",
      PlotDatas$Scenario[1] == "S11" ~ "Scenario 1 - Part 1: Model Selection via Information Criterion",
      PlotDatas$Scenario[1] == "S12" ~ "Scenario 1 - Part 2: Model then Variable Selection",
      PlotDatas$Scenario[1] == "S21" ~ "Scenario 2 - Part 1: Backward Variable Selection",
      PlotDatas$Scenario[1] == "S22" ~ "Scenario 2 - Part 2: Variable then Model Selection",
      TRUE ~ "Other"
    )
    
    Cor <- ifelse(PlotDatas$Correlation[1] == "Independent", "Independent", "Mod. Correlation")
    Cens <- ifelse(PlotDatas$Censoring[1] == "LowCens", "Low Censoring", "High Censoring")
    GenDist <- dplyr::case_when(
      PlotDatas$Gen.Distribution[1] == "weibull"     ~ "Weibull",
      PlotDatas$Gen.Distribution[1] == "exponential" ~ "Exponential",
      PlotDatas$Gen.Distribution[1] == "lognormal"   ~ "LogNormal",
      PlotDatas$Gen.Distribution[1] == "cox"         ~ "Gompertz",
      TRUE ~ PlotDatas$Gen.Distribution[1]
    )
    
    p <- ggplot(PlotDatas, aes(x = Regressors, y = Distribution, fill = Count)) +
      geom_tile(color = "white", linewidth = 0.75) +
      geom_text(aes(label = Count), color = "darkorange", size = text_size) +
      scale_fill_continuous_sequential(palette = "Dark Mint", limits = fill_range) +
      scale_x_discrete(limits = all_models) +
      scale_y_discrete(limits = all_distributions) +  # Ensure consistent y-axis across all plots
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio) +
      theme_classic(base_size = 8) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 8),
        plot.margin = margin(2, 2, 2, 2),  # Reducing plot margin to decrease white space
        legend.text = element_text(size = 5),     # Smaller text
        legend.title = element_text(size = 5),    # Smaller title
        legend.key.size = unit(1.25, "lines")     # Smaller key size
      )
    
    # Clean X-axis labels for all but bottom row
    if (i < 13) {
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    
    # Clean Y-axis labels for all but leftmost column
    if (!(i %in% c(1, 5, 9, 13))) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    
    if(i < 5) {
      p <- p + labs(title = paste(GenDist))
    }
    
    heatmap_plots[[i]] <- p
  }
  
  # ---- Create external row labels (left side) ----
  row_labels <- lapply(c(1, 5, 9, 13), function(i) {
    # Get the first plot in the row
    PlotDatas <- PlottingData[[i]]
    
    # Extract info for label
    Cor <- ifelse(PlotDatas$Correlation[1] == "Independent", "Independent", "Mod. Corr.")
    Cens <- ifelse(PlotDatas$Censoring[1] == "LowCens", "Low Cens.", "High Cens.")
    
    ggplot() + 
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste(Cens, "\n", Cor), 
               angle = 90, size = 3)
  })
  
  # Repeat row labels to match each row
  row_labels_full <- rep(row_labels, each = 1)
  
  # ---- Create external column labels (top) ----
  col_labels <- lapply(1:4, function(i) {
    ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = -5, 
               label = "", size = 0.5)  # 👈 Larger text
  })
  
  # ---- Group heatmaps by rows ----
  heatmap_rows <- split(heatmap_plots, ceiling(seq_along(heatmap_plots) / 4))
  
  # ---- Combine each row with its label on the left ----
  row_with_labels <- mapply(function(label, plots) {
    wrap_plots(label, wrap_plots(plots, nrow = 1), ncol = 2, widths = c(0.05, 1))  # 👈 Narrower spacing
  }, row_labels_full, heatmap_rows, SIMPLIFY = FALSE)
  
  # ---- Combine all rows vertically ----
  combined_body <- wrap_plots(row_with_labels, ncol = 1)
  
  # ---- Combine column labels with plot ----
  top_label_row <- wrap_plots(plot_spacer(), wrap_plots(col_labels, nrow = 1), 
                              ncol = 2, widths = c(-0.25, 0))  # 👈 Narrower spacing
  
  # ---- Combine everything ----
  final_plot <- wrap_plots(top_label_row, combined_body, ncol = 1, 
                           heights = c(0, 1), widths = c(0.05, 0.95)) +  # 👈 Tighter top spacing
    plot_layout(guides = "collect") & 
    plot_annotation(
      title = tittle,
      caption = Scen,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0, face = "italic")
        # legend.position = "none",                 # Top position for combined legend
        #legend.justification = c(1, 1)           # Align to the right within the top
      )
    )
  
  return(final_plot)
}
HeatmapPlot.Fun(Scenario11SubPlotList2, ratio = 0.10, text_size = 3,
                tittle = "Scenario 1 - Part 1: Model Selection Frequencies by Method and Simulation Setting")
HeatmapPlot.Fun(Scenario12SubPlotList2, ratio = 0.75, text_size = 2,
                tittle = "Scenario 1 - Part 2: Model Selection Frequencies by Method and Simulation Setting")
HeatmapPlot.Fun(Scenario21SubPlotList2, ratio = 0.75, text_size = 2,
                tittle = "Scenario 2 - Part 1: Model Selection Frequencies by Method and Simulation Setting")
HeatmapPlot.Fun(Scenario22SubPlotList2, ratio = 0.75, text_size = 2,
                tittle = "Scenario 2 - Part 2: Model Selection Frequencies by Method and Simulation Setting")

### Table Objects (Each scenario after marginalizing over regressor variable)
## Scenario 00 Table Set-Up
Scenario00_tbl <- bind_rows(Scenario00SubPlotList)

# Remove averaging and scenario column 
Scenario00_tbl[c("Averaging", "Scenario")] <- NULL

# Remove columns that end with lower, upper, true, count
Scenario00_tbl <- select(Scenario00_tbl, -c(matches("Lower"), matches("Upper"), matches("True"), matches("Count")))

# To rearrange fitted distribution, correlation, censoring, generating distribution to match the latex table creation 
Reorder_Tbl <- function(tbl){
  # Convert the columns to factors with the specified order
  tbl$Gen.Distribution <- factor(tbl$Gen.Distribution, 
                                            levels = c("weibull", "exponential", "lognormal", "cox"))
  
  tbl$Censoring <- factor(tbl$Censoring, 
                                     levels = c("LowCens", "HighCens"))
  
  tbl$Correlation <- factor(tbl$Correlation, 
                                       levels = c("Independent", "ModerateCorrelation"))
  
  tbl$Distribution <- factor(tbl$Distribution, 
                                        levels = c("AFT1_Weibull", "AFT2_Exponential", 
                                                   "AFT3_Lognormal", "Cox_Regression"))
  
  # Sort the data frame based on the specified order
  tbl <- tbl[order(tbl$Gen.Distribution, tbl$Censoring, tbl$Correlation, tbl$Distribution), ]
  # Return
  return(tbl)
}

# Do so
Scenario00_tbl <- Reorder_Tbl(Scenario00_tbl)

library(xtable)
#####
# xtable for Scenario00_tbl
#####
df <- data.frame(
  Distribution=rep(c(""), 64),
  Censoring=rep(c(""), 64),
  Correlation=rep(c(""), 64),
  Gen.Distribution=rep(c("Weibull", "Exponential", "Lognormal", "Gompertz"), 16),
  Beta1=round(Scenario00_tbl$Beta1, 3),
  Beta2=round(Scenario00_tbl$Beta2, 3),
  Beta3=round(Scenario00_tbl$Beta3, 3),
  Beta4=round(Scenario00_tbl$Beta4, 3),
  Beta_1_Bias=round(Scenario00_tbl$Beta_1_Bias, 3),
  Beta_2_Bias=round(Scenario00_tbl$Beta_2_Bias, 3),
  Beta_3_Bias=round(Scenario00_tbl$Beta_3_Bias, 3),
  Beta_4_Bias=round(Scenario00_tbl$Beta_4_Bias, 3),
  Beta1_Captured=round(Scenario00_tbl$Beta1_Captured, 3),
  Beta2_Captured=round(Scenario00_tbl$Beta2_Captured, 3),
  Beta3_Captured=round(Scenario00_tbl$Beta3_Captured, 3),
  Beta4_Captured=round(Scenario00_tbl$Beta4_Captured, 3)
)

# Formatting multi-row cells
# Generating Distribution Rows
row.lengths.1 <- rep(16, 4)
df[[1]][1 + (0:3) * 16] <- 
  paste0("\\midrule\\multirow{", row.lengths.1, "}{*}{\\textbf{", c("Weib", "Exp", "Lognorm", "Gomp"), "}}")

# Censoring Level Rows
row.lengths.2 <- rep(8, 8)
df[[2]][c(1, 9, 17, 25, 33, 41, 49, 57)] <- 
  paste0("\\multirow{", row.lengths.2, "}{*}{\\textbf{", rep(c("High", "Low"), 4), "}}")

# Correlation Level Rows
row.lengths.3 <- rep(4, 16)
df[[3]][c(1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61)] <- 
  paste0("\\multirow{", row.lengths.3, "}{*}{\\textbf{", rep(c("Ind", "Cor"), 8), "}}")


# Caption and xtable printing
strCaption <- paste0("\\textbf{Scenario 0 - Part 0:} Average Estimated Coefficients for Full Fitted Models Under Different Simulation Settings")

print(xtable(df, digits = c(0, 0, 0, 0, 0, rep(3, 12)), 
             align = "lllll|cccc|cccc|cccc",  # Decimnal align
             caption = strCaption, label = "tab:Table3"),
      size = "footnotesize",
      include.rownames = FALSE, 
      include.colnames = FALSE,
      caption.placement = "top", 
      hline.after = NULL, 
      floating = TRUE,
      sanitize.text.function = force, 
      add.to.row = list(pos = list(-1, nrow(df)),
                        command = c(paste("\\toprule\n",  
                                          "Generating & Censoring & \\multirow{2}{*}{Correlation} & Fitted & ",
                                          "\\multicolumn{4}{c|}{Average $\\beta$ Estimates} & ",
                                          "\\multicolumn{4}{c|}{Average $\\beta$ Bias} & ",
                                          "\\multicolumn{4}{c}{Average $\\beta$ Coverage} \\\\ \n", 
                                          "Distribution & Level & & Distribution & $\\beta_1$ & $\\beta_2$ & $\\beta_3$ & $\\beta_4$ & ",
                                          "$\\beta_1$ & $\\beta_2$ & $\\beta_3$ & $\\beta_4$ & ",
                                          "$\\beta_1$ & $\\beta_2$ & $\\beta_3$ & $\\beta_4$ \\\\ \n",
                                          "\\midrule\n"),
                                    paste("\\bottomrule\n")
                        )
      )
)

## Scenario 11 Table Set-Up
Scenario11_tbl <- bind_rows(Scenario11SubPlotList)

# Remove averaging and scenario columns and columns that end with lower, upper, true
Scenario11_tbl <- select(Scenario11_tbl, -c(matches("Lower"), matches("Upper"), 
                                            matches("True"),
                                            matches("Averaging"), matches("Scenario")))

# To rearrange settings to match the latex table creation 
Scenario11_tbl <- Reorder_Tbl(Scenario11_tbl)


### How to deal with dynamic plotting requirements...????


## Scenario 12 Table Set-Up
Scenario12_tbl <- bind_rows(Scenario12SubPlotList)

# Widening dataframe by Averaging column
# Separate the Conditional and Unconditional data
conditional_df <- subset(Scenario12_tbl, Averaging == "Conditional")
unconditional_df <- subset(Scenario12_tbl, Averaging == "Unconditional")

# Ensure that the columns related to beta values are the only ones that will be duplicated
# Get the columns that are common between both subsets (e.g., non-beta columns like Distribution, Scenario, etc.)
non_beta_columns <- setdiff(names(conditional_df), c("Averaging", "Beta1", "Beta2", "Beta3", "Beta4"))

# Create the wide table by combining the non-beta columns and adding the beta columns for both averaging types
Scenario12_tbl_wide <- cbind(
  conditional_df[, non_beta_columns], # Non-beta columns (kept only once)
  conditional_df[, grep("^Beta", names(conditional_df))], # Beta columns for Conditional
  unconditional_df[, grep("^Beta", names(unconditional_df))]  # Beta columns for Unconditional
)

# Rename the columns to clearly reflect the Averaging type
names(Scenario12_tbl_wide) <- c(
  non_beta_columns, # Non-beta columns (same for both Conditional and Unconditional)
  paste0("Conditional_", names(conditional_df)[grep("^Beta", names(conditional_df))]),  # Conditional beta columns
  paste0("Unconditional_", names(unconditional_df)[grep("^Beta", names(unconditional_df))]) # Unconditional beta columns
)

# Check the output
head(Scenario12_tbl_wide)


Scenario21_tbl <- bind_rows(Scenario21SubPlotList)
Scenario22_tbl <- bind_rows(Scenario22SubPlotList)
