### Sampling Distribution of Sample Variance
library(ggplot2)
Samp.Var <- function(nsim = 100000, n = 10){
  var.vec <- rep(0, nsim)
  for(i in 1:nsim){
    da.data <- rnorm(n, mean = 0, sd = 1)
    var.vec[i] <- var(da.data)
  }
  da.plot <- ggplot(data = as.data.frame(var.vec), aes(x = var.vec)) +
    stat_function(fun = function(x){dchisq(x*(n-1), df = n-1) * (n-1)},
                  color = "blue", size = 1.2) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.0125, fill = "#ce36f7", color = "black", alpha = 0.7) +
    xlim(min(var.vec), max(var.vec)) +
    labs(title = "Distribution of Sample Variances",
         x = "Sample Variance",
         y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12))
  print(da.plot)
}
Samp.Var(n = 5)