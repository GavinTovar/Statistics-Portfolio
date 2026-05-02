### Load Simulation Results
# load("./WNAR Paper/Simulation Results/CISubPlotData.RData")
# load("./WNAR Paper/Simulation Results/TabulatedResults.RData")

# Source functions and packages
source("./WNAR Paper/R Code/WNAR - Functions and Packages.R")

### Plotting
# Load Packages
library(tidyverse)
library(patchwork)
library(grid)
library(RColorBrewer)
library(colorspace)
# RColorBrewer::display.brewer.all()

# Some variables needed later...
nsim <- 30000

## Data manipulation
# To smoosh each table in the list into a single row (na.rm = TRUE)
Smooshed.Table.Results <- list()
for (i in 1:length(Tabulated.Results)) {
  Smooshing <- lapply(Tabulated.Results[[i]], FUN = function(x) {
    if (is.numeric(x)) {
      # Replace 0s with NAs (specifically for LASSO)
      x[x == 0] <- NA
      weighted.mean(x, w = Tabulated.Results[[i]][["Chosen.Freq"]], na.rm = TRUE)
    } else {
      x[1]
    }
  })
  
  Smooshed.Table.Results[[i]] <- bind_cols(Smooshing)
}

# Get rid of model and model frequency columns - Change from tibble to df
Smooshed.Table.Results <- as.data.frame(bind_rows(Smooshed.Table.Results)[, -c(1, 2)])

# Getting rid of correct and difference columns (Treat NAs as 0)
Smooshed.Table.Results <- Smooshed.Table.Results[, !grepl("(_C|_D)$", names(Smooshed.Table.Results))]

# Add a variable indicating which averaging was done
Smooshed.Table.Results[["Averaging"]] <- "NAsRemoved"
Smooshed.Table.Results

Smooshed.Table.Results2 <- list()

for(i in 1:length(Tabulated.Results)){
  Smooshing <- lapply(Tabulated.Results[[i]], FUN = function(x){
    if(is.numeric(x)){
      # Replace NAs with 0s
      x[is.na(x)] <- 0
      weighted.mean(x, w = Tabulated.Results[[i]][["Chosen.Freq"]])
    } else {
      x[1]
    }
  })
  Smooshed.Table.Results2[[i]] <- bind_cols(Smooshing)
}

# Get rid of model and model frequency columns - Change from tibble to df
Smooshed.Table.Results2 <- as.data.frame(bind_rows(Smooshed.Table.Results2)[, -c(1, 2)])

# Getting rid of correct and difference columns
Smooshed.Table.Results2 <- Smooshed.Table.Results2[, !grepl("(_C|_D)$", names(Smooshed.Table.Results2))]

# Add a variable indicating which averaging was done
Smooshed.Table.Results2[["Averaging"]] <- "NAs0"
Smooshed.Table.Results2

# Reshape Data for plotting and combine two dataframes into one
SmooshyTable <- bind_rows(Reshape.Fun(Smooshed.Table.Results), Reshape.Fun(Smooshed.Table.Results2))

# Combine model selection type and information criteria
SmooshyTable[["ModelSelectionType"]] <- ifelse(is.na(SmooshyTable[["IC"]]), 
                                               SmooshyTable[["MS.Type"]],
                                               paste(SmooshyTable[["MS.Type"]], SmooshyTable[["IC"]], sep = "_"))

### Subset data into different unique simulation settings
SubPlotDataList <- SmooshyTable %>%
  group_by(SampleSize, Error, Correlation, BetaSetting) %>%
  group_split()

# Fix bias column
SubPlotDataList <- lapply(SubPlotDataList, FUN = function(x) {
  x[["Bias"]] <- x[["Estimate"]] - x[["True"]]
  return(x)
})

### Line Plot Action
individual_plots <- list()
bias_range <- range(unlist(lapply(SubPlotDataList, function(df) df[["Bias"]])), na.rm = TRUE)

for (i in 1:(length(SubPlotDataList) - 1)) {
  PlotDatas <- SubPlotDataList[[i]]
  
  # Remove intercept coefficient rows
  PlotDatas <- subset(PlotDatas, Beta != "Beta0")
  
  # Naming stuff
  if(PlotDatas[["SampleSize"]][1] == "LowSS"){
    SS <- "Low Sample Size"
  } else {
    SS <- "High Sample Size"
  }
  
  if(PlotDatas[["Error"]][1] == "LowError"){
    Er <- "Low Error"
  } else {
    Er <- "High Error"
  }
  
  if(PlotDatas[["Correlation"]][1] == "Independent"){
    Cor <- "Independent"
  } else {
    Cor <- "Moderate Correlation"
  }
  
  if(PlotDatas[["BetaSetting"]][1] == "PatternA"){
    BetaSet <- "Coefficient Pattern A"
  } else {
    BetaSet <- "Coefficient Pattern B"
  }
  
  true_levels <- sort(unique(PlotDatas$True))
  custom_label <- function(x) {
    ifelse(x == 0, "0", formatC(x, format = "f", digits = 3, drop0trailing = TRUE))
  }
  
  individual_plots[[i]] <- ggplot(PlotDatas) +
    coord_cartesian(ylim = bias_range) +
    geom_line(aes(y = Bias, x = True, color = ModelSelectionType, linetype = Averaging)) +
    geom_point(aes(y = Bias, x = True, color = ModelSelectionType, shape = Averaging)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = true_levels,
                       labels = custom_label) +
    # scale_color_paletteer_d("ggsci::flattastic_flatui") +
    labs(y = paste(SS, "\n", Er, sep = ""),
         x = "Data Generating Coefficient Value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = margin(0, 0, 0, 0),
          axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 10))
  
  # Remove Y-axis labels for all but LHS plots
  if(!(i %in% c(1, 5, 9, 13))){
    individual_plots[[i]] <- individual_plots[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(y = "")
  }
  
  # Remove X-axis labels for all but bottom row
  if(i < 13){
    individual_plots[[i]] <- individual_plots[[i]] +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = "")
  }
  
  # Add subtitle to top row
  if(i < 5){
    individual_plots[[i]] <- individual_plots[[i]] +
      labs(subtitle = paste(Cor, "\n", BetaSet, sep = ""))
  }
}

# Combine all plots
combined_plot <- wrap_plots(individual_plots) + 
  plot_layout(ncol = 4, guides = "collect") &
  plot_annotation(
    title = "Post-Model Selection Bias (Coefficient Estimates - True Betas)",
    caption = "",
    theme = theme(
      plot.margin = margin(2, 2, 0.5, 2), # Top, right, bottom, left
      plot.caption = element_text(size = 10, face = "italic", hjust = 0, vjust = 3),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom"
    )
  )

print(combined_plot)

#### Heatmap Action
# Get rid of columns with coefficient values
Tabulated.Results.Trimmed <- lapply(Tabulated.Results, FUN = function(x) select(x, -matches("^Beta\\d+")))

# Add counts column
Tabulated.Results.Trimmed <- lapply(Tabulated.Results.Trimmed, FUN = function(x) {
  x[["Counts"]] <- x[["Chosen.Freq"]] * nsim 
  return(x)
})

# Stack all the dataframes
FrequencyFrame <- bind_rows(Tabulated.Results.Trimmed)

# Get rid of rows with Ridge
FrequencyFrame <- subset(FrequencyFrame, MS.Type != "Ridge")

# Combing Model selection type with information criteria
FrequencyFrame[["ModelSelectionType"]] <- ifelse(is.na(FrequencyFrame[["IC"]]), 
                                                 FrequencyFrame[["MS.Type"]],
                                                 paste(FrequencyFrame[["MS.Type"]], FrequencyFrame[["IC"]], sep = "_"))


# Optional: consistent model levels across all heatmaps
all_models <- sort(unique(FrequencyFrame$Model))

# Optional: consistent fill scale
fill_range <- range(FrequencyFrame$Counts, na.rm = TRUE)

# Split data by simulation setting
SubPlotDataList2 <- FrequencyFrame %>%
  group_by(SampleSize, Error, Correlation, BetaSetting) %>%
  group_split()

# ---- Generate heatmaps ----
heatmap_plots <- list()

for (i in seq_along(SubPlotDataList2)) {
  PlotDatas <- SubPlotDataList2[[i]]
  
  # Naming stuff
  if(PlotDatas[["SampleSize"]][1] == "LowSS"){
    SS <- "Low Sample Size"
  } else {
    SS <- "High Sample Size"
  }
  
  if(PlotDatas[["Error"]][1] == "LowError"){
    Er <- "Low Error"
  } else {
    Er <- "High Error"
  }
  
  if(PlotDatas[["Correlation"]][1] == "Independent"){
    Cor <- "Independent"
  } else {
    Cor <- "Moderate Correlation"
  }
  
  if(PlotDatas[["BetaSetting"]][1] == "PatternA"){
    BetaSet <- "Coefficient Pattern A"
  } else {
    BetaSet <- "Coefficient Pattern B"
  }
  
  p <- ggplot(PlotDatas, aes(x = Model, y = ModelSelectionType, fill = Counts)) +
    geom_tile(color = "white", linewidth = 0.75) +
    scale_fill_continuous_sequential(palette = "SunsetDark", limits = fill_range) +
    scale_x_discrete(limits = all_models) +
    labs(x = NULL, y = NULL) +
    coord_fixed(2) +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      panel.border = element_rect(color = "black", fill = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Remove Y-axis labels for non-left plots
  if (!(i %in% c(1, 5, 9, 13))) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  heatmap_plots[[i]] <- p
}

# ---- Pad to 4x4 if needed (add empty spaces) ----
while (length(heatmap_plots) < 16) {
  heatmap_plots[[length(heatmap_plots) + 1]] <- plot_spacer()
}

# ---- Create external row labels (left side) ----
row_labels <- lapply(c(1, 5, 9, 13), function(i) {
  PlotDatas <- SubPlotDataList2[[i]]
  if(PlotDatas[["SampleSize"]][1] == "LowSS"){
    SS <- "Low Sample Size"
  } else {
    SS <- "High Sample Size"
  }
  
  if(PlotDatas[["Error"]][1] == "LowError"){
    Er <- "Low Error"
  } else {
    Er <- "High Error"
  }
  
  ggplot() + 
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, 
             label = paste(SS, "\n", Er, sep = ""), 
             angle = 90, size = 4.5)  # 👈 Larger text
})

# Repeat row labels to match each row
row_labels_full <- rep(row_labels, each = 1)

# ---- Create external column labels (top) ----
col_labels <- lapply(1:4, function(i) {
  PlotDatas <- SubPlotDataList2[[i]]
  if(PlotDatas[["Correlation"]][1] == "Independent"){
    Cor <- "Independent"
  } else {
    Cor <- "Moderate Correlation"
  }
  
  if(PlotDatas[["BetaSetting"]][1] == "PatternA"){
    BetaSet <- "Coefficient Pattern A"
  } else {
    BetaSet <- "Coefficient Pattern B"
  }
  
  ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, 
             label = paste(Cor, "\n", BetaSet, sep = ""), size = 4.5)  # 👈 Larger text
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
                            ncol = 2, widths = c(0.075, 1))  # 👈 Narrower spacing

# ---- Combine everything ----
final_plot <- wrap_plots(top_label_row, combined_body, ncol = 1, heights = c(0.05, 1)) +  # 👈 Tighter top spacing
  plot_annotation(
    title = "Model Selection Frequencies by Method and Simulation Setting",
    caption = "",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0, face = "italic")
    )
  )

# ---- Show the plot ----
print(final_plot)

### Coverage Plot
# Adding True beta columns
for (i in 1:length(CISubPlotDataList)) {
  if(CISubPlotDataList[[i]][["BetaSetting"]][1] == "PatternA"){
    CISubPlotDataList[[i]][["True"]] <- rep(c(0, c(0.125, 0.25, 0.50, 1, 2)), 12)
  } else {
    CISubPlotDataList[[i]][["True"]] <- rep(c(0, c(-1.5, -0.75, 0, 0.75, 1.5)), 12)
  }
}

# Combining model selection with information criteria
CISubPlotDataList <- lapply(CISubPlotDataList, FUN = function(x){
  x[["ModelSelectionType"]] <- ifelse(is.na(x[["IC"]]), 
                                                 x[["MS.Type"]],
                                                 paste(x[["MS.Type"]], x[["IC"]], sep = "_"))
  return(x)
})

# Interleave the two halves without combining them
CISubPlotData <- CISubPlotDataList
n <- length(CISubPlotDataList) / 2
CISubPlotData[seq(1, n * 2, by = 2)] <- CISubPlotDataList[1:n]           # A's into odd positions
CISubPlotData[seq(2, n * 2, by = 2)] <- CISubPlotDataList[(n + 1):(n*2)] # B's into even positions

# Starting plotting
individual_plots <- list()

for (i in 1:length(CISubPlotData)) {
  PlotDatas <- CISubPlotData[[i]]
  
  # Remove intercept coefficient rows
  PlotDatas <- subset(PlotDatas, Beta != "Beta0")
  
  # Naming stuff
  if(PlotDatas[["SampleSize"]][1] == "LowSS"){
    SS <- "Low Sample Size"
  } else {
    SS <- "High Sample Size"
  }
  
  if(PlotDatas[["Error"]][1] == "LowError"){
    Er <- "Low Error"
  } else {
    Er <- "High Error"
  }
  
  if(PlotDatas[["Correlation"]][1] == "Independent"){
    Cor <- "Independent"
  } else {
    Cor <- "Moderate Correlation"
  }
  
  if(PlotDatas[["BetaSetting"]][1] == "PatternA"){
    BetaSet <- "Coefficient Pattern A"
  } else {
    BetaSet <- "Coefficient Pattern B"
  }
  
  true_levels <- sort(unique(PlotDatas$True))
  custom_label <- function(x) {
    ifelse(x == 0, "0", formatC(x, format = "f", digits = 3, drop0trailing = TRUE))
  }
  
  individual_plots[[i]] <- ggplot(PlotDatas) +
    geom_line(aes(y = Coverage, x = True, color = ModelSelectionType, linetype = Averaging)) +
    geom_point(aes(y = Coverage, x = True, color = ModelSelectionType, shape = Averaging)) +
    geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = true_levels,
                       labels = custom_label) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Set3", aesthetics = "color") +
    labs(y = paste(SS, "\n", Er, sep = ""),
         x = "Data Generating Coefficient Value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = margin(0, 0, 0, 0),
          axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 10))
  
  # Remove Y-axis labels for all but LHS plots
  if(!(i %in% c(1, 5, 9, 13))){
    individual_plots[[i]] <- individual_plots[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(y = "")
  }
  
  # Remove X-axis labels for all but bottom row
  if(i < 13){
    individual_plots[[i]] <- individual_plots[[i]] +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = "")
  }
  
  # Add subtitle to top row
  if(i < 5){
    individual_plots[[i]] <- individual_plots[[i]] +
      labs(subtitle = paste(Cor, "\n", BetaSet, sep = ""))
  }
}

# Combine all plots
combined_plot <- wrap_plots(individual_plots) + 
  plot_layout(nrow = 4, guides = "collect") &
  plot_annotation(
    title = "Coefficient Coverage",
    caption = "",
    theme = theme(
      plot.margin = margin(2, 2, 0.5, 2), # Top, right, bottom, left
      plot.caption = element_text(size = 10, face = "italic", hjust = 0, vjust = 3),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom"
    )
  )

print(combined_plot)

### PCS - Like Table
library(xtable)

# Assuming tbldata is already available as shown
tbldata <- select(Tabulated.Results[[1]], -c(matches("_B"), matches("_C"), matches("_D"),
                                             matches("SampleSize"),
                                             matches("Correlation"), matches("Error"),
                                             matches("MS.Type"), matches("IC"), 
                                             matches("BetaSetting")))
# Rearrange to match plotting
# Reorder the columns in the desired sequence
tbldata <- tbldata[, c("Model", "Chosen.Freq", 
                       "Beta0", "Beta0_E", "Beta1", "Beta1_E",
                       "Beta2", "Beta2_E", "Beta3", "Beta3_E",
                       "Beta4", "Beta4_E", "Beta5", "Beta5_E")]

# Number of columns in tbldata (20) plus 1 for row names
n_cols <- ncol(tbldata) + 1

# Construct the alignment string: "l" for row names, "l" for first two columns, "c" for remaining
align_str <- paste(c("l", "l", "l", rep("c", n_cols - 3)), collapse = "")

# align_str <- paste(c("l", 
#                      "|>{\\columncolor{lightgray}}l", 
#                      ">{\\columncolor{lightgray}}l", 
#                      rep(">{\\columncolor{lightpurple}}c", 6), 
#                      rep(">{\\columncolor{lightblue}}c", 6)), collapse = "")

# Caption for the table
strCaption <- paste0("\\textbf{Caption Placeholder}")

# Print the xtable with customized alignment and header
print(xtable(tbldata, digits = c(0, 0, 3, rep(3, n_cols - 3)), 
             align = align_str, 
             caption = strCaption, label = "tab:Table3"),
      size = "footnotesize",
      include.rownames = FALSE, 
      include.colnames = FALSE,
      caption.placement = "top", 
      hline.after = NULL, 
      floating = TRUE,
      sanitize.text.function = force, 
      add.to.row = list(pos = list(-1, nrow(tbldata)),
                        command = c(paste("\\toprule\n",  
                                          "\\multicolumn{2}{c|}{Model Info} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{0} = 0$} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{1} = -1.5$} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{2} = -0.75$} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{3} = 0$} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{4} = 0.75$} & ",
                                          "\\multicolumn{2}{c|}{$\\beta_{5} = 1.5$} \\\\ \n",
                                          "Model & Chosen.Freq & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ & ",
                                          "$\\hat{\\beta}_{PS}$ & $\\hat{\\beta}^{*}$ \\\\ \n",
                                          "\\midrule\n"),
                                    paste("\\bottomrule\n")
                        )
      )
)
