### (ST 599) Homework 2
setwd("C:/Users/Ghcto/OneDrive/Desktop/School/Year 2/Winter 2025/ST 599/Homeworks")
# setwd("C:/Users/Ghcto/OneDrive/Desktop/School/Year 2/Winter 2025/ST 599/Homworks/HW2")

## Question 1 - Part 1
# Set seed for reproducibility
set.seed(529461)

# Define parameters
I <- 10
n <- 10
N <- I * n

sigma2 <- 0.5
tau2 <- 1
mu <- 4
phi0 <- -2
phi1 <- 0.25

# Generate group-level effects
alpha <- rnorm(I, mean = mu, sd = sqrt(tau2))

# Generate responses and missingness indicators
y <- numeric(N)
m <- numeric(N)
group <- rep(1:I, each = n)  # Group indices

for (k in 1:N) {
  y[k] <- rnorm(1, mean = alpha[group[k]], sd = sqrt(sigma2))
  pi_k <- 1 / (1 + exp(-(phi0 + phi1 * alpha[group[k]])))
  m[k] <- rbinom(1, 1, pi_k)
}

# Extract observed y values (where m == 0)
y_obs <- y[m == 0]
g_obs <- group[m == 0] # Group indices for observed values
N0 <- length(y_obs)  # Number of observed y values

# Prepare data list for Stan
StanData <- list(
  I = I,
  N = N,
  N0 = N0,
  group = group,
  g_obs = g_obs,
  y_obs = y_obs,
  m = m
)
StanData

## Question 1 - Part 2
library(knitr)
library(cmdstanr)

# Compile and fit the model
Q1P2.Model <- cmdstan_model("Question1Part2.stan")
Q1P2.Fit <- Q1P2.Model[["sample"]](data = StanData, seed = 529461, 
                                   chains = 4, parallel_chains = 4)
Q1P2.Fit[["summary"]]()

# LaTeX Table of results
kable(Q1P2.Fit[["summary"]](), format = "latex", booktabs = TRUE,
      caption = "Stan Model Fit Results")

## Question 1 - Part 3
# Compile and fit the model
Q1P3.Model <- cmdstan_model("Question1Part3.stan")
Q1P3.Fit <- Q1P3.Model[["sample"]](data = StanData, seed = 529461, 
                                   chains = 4, parallel_chains = 4)
Q1P3.Fit[["summary"]]()

# LaTeX Table of results
kable(Q1P3.Fit[["summary"]](), format = "latex", booktabs = TRUE,
      caption = "Stan Model Fit Results")

## Question 2
library(Surrogate)
data("Schizo_PANSS")
# Schizo_PANSS <- read.csv("Schizo_PANSS.csv")

hw_data <- Schizo_PANSS[,c("Id","Treat","Week1","Week2","Week4","Week6","Week8")]
hw_data_sub <- hw_data |>
  subset(Id %in% 1:200)

long_case <- stats::reshape(
  hw_data_sub,
  direction = "long",
  varying = 3:7,
  sep = ""
)[,-5]
names(long_case) <- c("Id","Treat","time","panss")

# Putting it into a list

# Create a missingness column for the response
# Create a vector of zeros of length 1000
long_case[["na_positions"]] <- rep(0, nrow(long_case))

# Get the positions of NA values
NA.Pos <- which(is.na(long_case[["panss"]]))

# Set the corresponding positions in the vector to 1
long_case[["na_positions"]][NA.Pos] <- 1

# Group indices
Group <- rep(1:length(unique(long_case[["Id"]])), 
             times = 5)

# Get indices of observed and missing values
obs_idx <- which(!is.na(long_case[["panss"]]))
miss_idx <- which(is.na(long_case[["panss"]]))

# Prepare data list for Stan
long_data <- list(
  I = length(unique(long_case[["Id"]])),    # Number of groups
  N = nrow(long_case),                      # Total number of observations
  N0 = length(obs_idx),                     # Number of observed y values
  group = Group,                            # Group indices
  g_obs = Group[obs_idx],                   # Observed Group indices
  y_obs = long_case[["panss"]][obs_idx],    # Observed y values
  m = long_case[["na_positions"]]           # Missingness indicator (0 = observed, 1 = missing)
)
long_data

# Compile and fit the model
Q2.Model <- cmdstan_model("Question2.stan")
Q2.Fit <- Q2.Model[["sample"]](data = long_data, seed = 529461, 
                               chains = 4, parallel_chains = 4)
Q2.Fit[["summary"]]()

# LaTeX Table of results
kable(Q2.Fit[["summary"]](), format = "latex", booktabs = TRUE,
      caption = "Stan Model Fit Results")

## Question 3 - Part 1
# Data structure from long_case
N <- nrow(long_case)
idx_i <- long_case$Id |> as.factor() |> as.integer()
I <- max(idx_i)
treat <- long_case$Treat
tilde_t <- scale(long_case$time)

# Define parameters
sigma2 <- 0.5
tau2_a <- 1
mu_a <- 4
phi0 <- -2
phi1 <- 0.25
tau2_b <- 8^2
mu_b <- -5
betatreat <- -0.5

# Making up some values for true phis
phi2 <- 0.5
phi3 <- 0.5

# Generate group-level effects
alpha <- rnorm(I, mean = mu_a, sd = sqrt(tau2_a))
beta <- rnorm(I, mean = mu_b, sd = sqrt(tau2_b))

# Generate responses and missingness indicators
y <- numeric(N)
m <- numeric(N)
group <- rep(1:I, each = N/I)  # Group indices

for (k in 1:N) {
  y[k] <- rnorm(1, mean = alpha[group[k]] + beta[group[k]] * tilde_t[k] + betatreat * treat[group[k]], sd = sqrt(sigma2))
  pi_k <- 1 / (1 + exp(-(phi0 + phi1 * alpha[group[k]] + phi2 * beta[group[k]] + phi3 * treat[group[k]])))
  m[k] <- rbinom(1, 1, pi_k)
}

# Extract observed y values (where m == 0)
y_obs <- y[m == 0]     # Indices of observed y values
g_obs <- group[m == 0] # Group indices for observed values
N0 <- length(y_obs)    # Number of observed y values

# Prepare data list for Stan
StanData3 <- list(
  I = I,
  N = N,
  N0 = N0,
  tilde_t = as.vector(tilde_t),
  tilde_t_obs = as.vector(tilde_t[m == 0]),
  group = group,
  g_obs = g_obs,
  y_obs = y_obs,
  m = m,
  treat = treat,
  treat_obs = treat[m == 0]
)
StanData3

# Compile and fit the model
Q3P1.Model <- cmdstan_model("Question3Part1.stan")
Q3P1.Fit <- Q3P1.Model[["sample"]](data = StanData3, seed = 529461, 
                                   chains = 4, parallel_chains = 4)
Q3P1.Fit[["summary"]]()

# LaTeX Table of results
kable(Q3P1.Fit[["summary"]](), format = "latex", booktabs = TRUE,
      caption = "Stan Model Fit Results")

# Prepare data list for Stan
Long_data <- list(
  I = length(unique(long_case[["Id"]])),    # Number of groups
  N = nrow(long_case),                      # Total number of observations
  N0 = length(obs_idx),                     # Number of observed y values
  tilde_t = as.vector(tilde_t),             # Scaled times
  group = group,                            # Group indices
  g_obs = group[obs_idx],                   # Observed Group indices
  y_obs = long_case[["panss"]][obs_idx],    # Observed y values
  m = long_case[["na_positions"]],          # Missingness indicator (0 = observed, 1 = missing)
  treat = treat                             # Treatment indicator
)
Long_data

# Compile and fit the model
Q3P2.Model <- cmdstan_model("Question3Part1.stan")
Q3P2.Fit <- Q3P2.Model[["sample"]](data = Long_data, seed = 529461, 
                                   chains = 4, parallel_chains = 4)
Q3P2.Fit[["summary"]]()