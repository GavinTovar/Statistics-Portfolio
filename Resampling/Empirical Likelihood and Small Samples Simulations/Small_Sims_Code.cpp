// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm> // for std::sort
using namespace Rcpp;

// Some notes

// Why a clone() is sometimes needed here
// Rcpp types (e.g., NumericVector) are "smart" wrappers around R's SEXP objects.
// These SEXPs point to memory owned by R.
// When you do: 
// NumericVector sorted = x;
// you’re just making another reference to the same R object — no deep copy happens


// [[Rcpp::export]]
NumericVector z_interval(NumericVector x, double parm, double alpha = 0.05) {
  int n = x.size();
  double mean_x = mean(x);
  double sd_x = sd(x);
  
  double z = R::qnorm(1 - alpha / 2, 0, 1, 1, 0); // lower.tail=TRUE, log=FALSE
  double lower = mean_x - z * sd_x / sqrt(n);
  double upper = mean_x + z * sd_x / sqrt(n);
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat z_interval_grid(NumericVector x, double parm, 
                          double alpha_start = 0.0005, double alpha_end = 0.20,
                          int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  double mean_x = mean(x);
  double sd_x = sd(x);
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    double z = R::qnorm(1 - alpha / 2, 0, 1, 1, 0);
    double lower = mean_x - z * sd_x / sqrt(n);
    double upper = mean_x + z * sd_x / sqrt(n);
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// [[Rcpp::export]]
NumericVector t_interval(NumericVector x, double parm, double alpha = 0.05) {
  int n = x.size();
  double mean_x = mean(x);
  double sd_x = sd(x);
  int df = n - 1;
  
  double t = R::qt(1 - alpha / 2, df, 1, 0); // lower.tail=TRUE, log=FALSE
  double lower = mean_x - t * sd_x / sqrt(n);
  double upper = mean_x + t * sd_x / sqrt(n);
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat t_interval_grid(NumericVector x, double parm, 
                          double alpha_start = 0.0005, double alpha_end = 0.20,
                          int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  double mean_x = mean(x);
  double sd_x = sd(x);
  int df = n - 1;
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    double t = R::qt(1 - alpha / 2, df, 1, 0);
    double lower = mean_x - t * sd_x / sqrt(n);
    double upper = mean_x + t * sd_x / sqrt(n);
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// Helper function to compute an empirical quantile of a numeric vector
// Arguments:
//   x - NumericVector (Rcpp) containing data values
//   p - quantile level between 0 and 1 (e.g., 0.5 for median)
// Returns: The empirical p-th quantile, computed via linear interpolation
double empiricalQuantile(NumericVector x, double p) {
  // Create a copy of the input vector so the original is not modified
  NumericVector sorted = clone(x);
  
  // Sort the copied vector in ascending order
  std::sort(sorted.begin(), sorted.end());
  
  // Get the number of elements in the sorted vector
  int n = sorted.size();
  
  // Compute the exact (floating-point) index position of the desired quantile
  // "n - 1" ensures the position is between the first (0) and last (n-1) index
  double pos = p * (n - 1);
  
  // Determine the lower index (integer) at or below the target position
  int lo = floor(pos);
  
  // Determine the upper index (integer) at or above the target position
  int hi = ceil(pos);
  
  // Fractional distance between lo and hi (used for interpolation)
  double h = pos - lo;
  
  // Perform linear interpolation between sorted[lo] and sorted[hi]
  // If h = 0, this is exactly sorted[lo]; if h = 1, exactly sorted[hi]
  return (1 - h) * sorted[lo] + h * sorted[hi];
}

// Function to compute percentile bootstrap confidence interval
// [[Rcpp::export]]
NumericVector perc_boot_interval(NumericVector x, double parm, int B = 100, double alpha = 0.05) {
  int n = x.size();
  NumericVector results(B);
  
  for(int b = 0; b < B; ++b) {
    double sum = 0.0;
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      sum += x[idx];
    }
    results[b] = sum / n;
  }
  
  NumericVector CI(4);
  CI[0] = empiricalQuantile(results, alpha / 2);
  CI[1] = empiricalQuantile(results, 1 - alpha / 2);
  
  CI[2] = (CI[0] < parm && parm < CI[1]) ? 1.0 : 0.0;
  CI[3] = CI[1] - CI[0];
  
  return CI;
}

// [[Rcpp::export]]
arma::mat perc_boot_interval_grid(NumericVector x, double parm, int B = 100,
                          double alpha_start = 0.0005, double alpha_end = 0.20,
                          int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  
  // Generate B bootstrap statistics
  NumericVector results(B);
  
  for(int b = 0; b < B; ++b) {
    double sum = 0.0;
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      sum += x[idx];
    }
    results[b] = sum / n;
  }
  
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    // Compute the lower and upper quantiles for the current alpha
    double lower = empiricalQuantile(results, alpha / 2);
    double upper = empiricalQuantile(results, 1 - alpha / 2);
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// [[Rcpp::export]]
NumericVector bca_boot_interval(NumericVector x, double parm, int B = 1000, double alpha = 0.05) {
  // Basic checks
  int n = x.size();
  if (n <= 1) {
    stop("x must have at least 2 observations for jackknife & bootstrap.");
  }
  if (B <= 0) {
    stop("B must be positive.");
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).");
  }
  
  // Observed statistic (theta_hat) - here we use the sample mean to match other functions
  double theta_hat = mean(x);

  // 1) Generate bootstrap replicates theta*_b (using simple nonparametric bootstrap)
  NumericVector theta_star(B);
  for (int b = 0; b < B; ++b) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n)); // uniform on [0,n), convert to index 0..n-1
      sum += x[idx];
    }
    theta_star[b] = sum / n;
  }
  
  // 2) Compute bias-correction factor z0:
  // proportion of bootstrap statistics less than the observed theta_hat
  int count_less = 0;
  for (int b = 0; b < B; ++b) if (theta_star[b] < theta_hat) ++count_less;
  double prop = double(count_less) / double(B);

  // clip prop to avoid exact 0 or 1 which give -Inf/+Inf when inverted
  double eps = 1e-12;
  if (prop < eps) prop = eps;
  if (prop > 1.0 - eps) prop = 1.0 - eps;

  double z0 = R::qnorm(prop, 0.0, 1.0, 1, 0); // inverse normal CDF

  // 3) Compute acceleration a using jackknife leave-one-out estimates
  // Efficiently compute leave-one-out means using total sum
  double total_sum = 0.0;
  for (int i = 0; i < n; ++i) total_sum += x[i];

  NumericVector theta_minus(n);
  for (int i = 0; i < n; ++i) {
    theta_minus[i] = (total_sum - x[i]) / double(n - 1); // leave-one-out mean
  }
  
  double theta_bar = mean(theta_minus);

  // compute numerator and denominator of acceleration formula
  double num = 0.0;
  double denom_part = 0.0;
  for (int i = 0; i < n; ++i) {
    double diff = theta_bar - theta_minus[i];
    num += diff * diff * diff;        // sum (theta_bar - theta_-i)^3
    denom_part += diff * diff;        // sum (theta_bar - theta_-i)^2
  }
  
  double a = 0.0;
  double denom = 6.0 * std::pow(denom_part, 1.5);
  if (denom_part <= 0.0 || denom == 0.0 || !R_FINITE(denom)) {
    // If denom_part is zero (no variability), acceleration is set to 0
    a = 0.0;
  } else {
    a = num / denom;
  }
  
  // 4) Compute adjusted quantiles alpha1 and alpha2
  // z_{alpha/2} and z_{1-alpha/2}
  double z_lower = R::qnorm(alpha / 2.0, 0.0, 1.0, 1, 0);
  double z_upper = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);
  
  // Helper lambda to compute adjusted percentile
  auto adjusted_alpha = [&](double z_alpha) {
    double denom_adj = 1.0 - a * (z0 + z_alpha);
    // Avoid divide-by-zero; if denom_adj is zero or near zero, push it away slightly
    if (std::abs(denom_adj) < 1e-12) {
      denom_adj = (denom_adj >= 0) ? 1e-12 : -1e-12;
    }
    double w = z0 + (z0 + z_alpha) / denom_adj;
    double res = R::pnorm(w, 0.0, 1.0, 1, 0); // Phi(w)
    // Clip to [0,1]
    if (res < 0.0) res = 0.0;
    if (res > 1.0) res = 1.0;
    return res;
  };
  
  double alpha1 = adjusted_alpha(z_lower);
  double alpha2 = adjusted_alpha(z_upper);
  
  // Clip alpha1/alpha2 to avoid extreme exact 0/1 (empiricalQuantile can't handle p exactly 0/1 well)
  double tiny = 1.0 / double(B + 1); // a small but reasonable step
  if (alpha1 < tiny) alpha1 = tiny;
  if (alpha2 < tiny) alpha2 = tiny;
  if (alpha1 > 1.0 - tiny) alpha1 = 1.0 - tiny;
  if (alpha2 > 1.0 - tiny) alpha2 = 1.0 - tiny;
  
  // Ensure alpha1 <= alpha2 (numerical safety)
  if (alpha1 > alpha2) std::swap(alpha1, alpha2);
  
  // 5) Get the bootstrap quantiles at adjusted percentiles
  double lower = empiricalQuantile(theta_star, alpha1);
  double upper = empiricalQuantile(theta_star, alpha2);
  
  // 6) capture indicator and length
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}


// [[Rcpp::export]]
arma::mat bca_boot_interval_grid(NumericVector x, double parm, int B = 100,
                                  double alpha_start = 0.0005, double alpha_end = 0.20,
                                  int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Basic checks
  int n = x.size();
  if (n <= 1) {
    stop("x must have at least 2 observations for jackknife & bootstrap.");
  }
  if (B <= 0) {
    stop("B must be positive.");
  }
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  
  
  // Observed statistic (theta_hat) - here we use the sample mean to match other functions
  double theta_hat = mean(x);
  
  // 1) Generate bootstrap replicates theta*_b (using simple nonparametric bootstrap)
  NumericVector theta_star(B);
  for (int b = 0; b < B; ++b) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n)); // uniform on [0,n), convert to index 0..n-1
      sum += x[idx];
    }
    theta_star[b] = sum / n;
  }
  
  // 2) Compute bias-correction factor z0:
  // proportion of bootstrap statistics less than the observed theta_hat
  int count_less = 0;
  for (int b = 0; b < B; ++b) if (theta_star[b] < theta_hat) ++count_less;
  double prop = double(count_less) / double(B);
  
  // clip prop to avoid exact 0 or 1 which give -Inf/+Inf when inverted
  double eps = 1e-12;
  if (prop < eps) prop = eps;
  if (prop > 1.0 - eps) prop = 1.0 - eps;
  
  double z0 = R::qnorm(prop, 0.0, 1.0, 1, 0); // inverse normal CDF
  
  // 3) Compute acceleration a using jackknife leave-one-out estimates
  // Efficiently compute leave-one-out means using total sum
  double total_sum = 0.0;
  for (int i = 0; i < n; ++i) total_sum += x[i];
  
  NumericVector theta_minus(n);
  for (int i = 0; i < n; ++i) {
    theta_minus[i] = (total_sum - x[i]) / double(n - 1); // leave-one-out mean
  }
  
  double theta_bar = mean(theta_minus);
  
  // compute numerator and denominator of acceleration formula
  double num = 0.0;
  double denom_part = 0.0;
  for (int i = 0; i < n; ++i) {
    double diff = theta_bar - theta_minus[i];
    num += diff * diff * diff;        // sum (theta_bar - theta_-i)^3
    denom_part += diff * diff;        // sum (theta_bar - theta_-i)^2
  }
  
  double a = 0.0;
  double denom = 6.0 * std::pow(denom_part, 1.5);
  if (denom_part <= 0.0 || denom == 0.0 || !R_FINITE(denom)) {
    // If denom_part is zero (no variability), acceleration is set to 0
    a = 0.0;
  } else {
    a = num / denom;
  }
  
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    // 4) Compute adjusted quantiles alpha1 and alpha2
    // z_{alpha/2} and z_{1-alpha/2}
    double z_lower = R::qnorm(alpha / 2.0, 0.0, 1.0, 1, 0);
    double z_upper = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);
    
    // Helper lambda to compute adjusted percentile
    auto adjusted_alpha = [&](double z_alpha) {
      double denom_adj = 1.0 - a * (z0 + z_alpha);
      // Avoid divide-by-zero; if denom_adj is zero or near zero, push it away slightly
      if (std::abs(denom_adj) < 1e-12) {
        denom_adj = (denom_adj >= 0) ? 1e-12 : -1e-12;
      }
      double w = z0 + (z0 + z_alpha) / denom_adj;
      double res = R::pnorm(w, 0.0, 1.0, 1, 0); // Phi(w)
      // Clip to [0,1]
      if (res < 0.0) res = 0.0;
      if (res > 1.0) res = 1.0;
      return res;
    };
    
    double alpha1 = adjusted_alpha(z_lower);
    double alpha2 = adjusted_alpha(z_upper);
    
    // Clip alpha1/alpha2 to avoid extreme exact 0/1 (empiricalQuantile can't handle p exactly 0/1 well)
    double tiny = 1.0 / double(B + 1); // a small but reasonable step
    if (alpha1 < tiny) alpha1 = tiny;
    if (alpha2 < tiny) alpha2 = tiny;
    if (alpha1 > 1.0 - tiny) alpha1 = 1.0 - tiny;
    if (alpha2 > 1.0 - tiny) alpha2 = 1.0 - tiny;
    
    // Ensure alpha1 <= alpha2 (numerical safety)
    if (alpha1 > alpha2) std::swap(alpha1, alpha2);
    
    // 5) Get the bootstrap quantiles at adjusted percentiles
    double lower = empiricalQuantile(theta_star, alpha1);
    double upper = empiricalQuantile(theta_star, alpha2);
    
    // 6) capture indicator and length
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// [[Rcpp::export]]
NumericVector bootstrap_t(NumericVector x, double parm, int B = 100, double alpha = 0.05) {
  int n = x.size();
  double theta_hat = mean(x);
  double se_hat = sd(x) / sqrt(n);
  
  NumericVector t_star(B);
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));  // bootstrap index
      x_star[i] = x[idx];
    }
    double theta_star = mean(x_star);
    double se_star = sd(x_star) / sqrt(n);
    t_star[b] = (theta_star - theta_hat) / se_star;
  }
  
  double lower_q = empiricalQuantile(t_star, alpha/2);
  double upper_q = empiricalQuantile(t_star, 1 - alpha/2);
  
  double lower = theta_hat - upper_q * se_hat;
  double upper = theta_hat - lower_q * se_hat;
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat bootstrap_t_grid(NumericVector x, double parm, int B = 100,
                          double alpha_start = 0.0005, double alpha_end = 0.20,
                          int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  double theta_hat = mean(x);
  double se_hat = sd(x) / sqrt(n);
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  
  // Bootstrap stats
  NumericVector t_star(B);
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      x_star[i] = x[idx];
    }
    double theta_star = mean(x_star);
    double se_star = sd(x_star) / sqrt(n);
    t_star[b] = (theta_star - theta_hat) / se_star;
  }
  
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    // Compute quantiles
    double lower_q = empiricalQuantile(t_star, alpha/2);
    double upper_q = empiricalQuantile(t_star, 1 - alpha/2);
    
    double lower = theta_hat - upper_q * se_hat;
    double upper = theta_hat - lower_q * se_hat;
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}


// Helper function to find the shortest interval containing (1 - alpha) fraction of values
// Assumes `sorted_x` is already sorted in ascending order
// Arguments:
//   sorted_x - sorted numeric values (Rcpp NumericVector)
//   alpha     - significance level (e.g., 0.05 for a 95% coverage interval)
// Returns: A std::pair {lower_bound, upper_bound} of the shortest such interval
std::pair<double, double> shortest_interval(NumericVector sorted_x, double alpha) {
  
  // Total number of observations
  int n = sorted_x.size();
  
  // Number of points that should be covered in the interval
  // ceil() is used in case (1-alpha)*n is not an integer
  int cover = (int)std::ceil((1.0 - alpha) * n);
  
  // Initialize minimum width to positive infinity
  double min_width = R_PosInf;
  
  // Index where the shortest interval starts (to be determined)
  int min_idx = 0;
  
  // Slide a window of length `cover` through the sorted vector
  for (int i = 0; i + cover - 1 < n; ++i) {
    
    // Width of the current interval:
    // last element in the window minus the first
    double width = sorted_x[i + cover - 1] - sorted_x[i];
    
    // If this interval is shorter than the best found so far, update the record
    if (width < min_width) {
      min_width = width;
      min_idx = i;
    }
  }
  
  // Return the bounds of the shortest interval found
  return std::make_pair(sorted_x[min_idx], 
                        sorted_x[min_idx + cover - 1]);
}

// [[Rcpp::export]]
NumericVector short_bootstrap_t(NumericVector x, double parm, int B = 100, double alpha = 0.05) {
  int n = x.size();
  double theta_hat = mean(x);
  double se_hat = sd(x) / sqrt(n);
  
  NumericVector t_star(B);
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      x_star[i] = x[idx];
    }
    double theta_star = mean(x_star);
    double se_star = sd(x_star) / sqrt(n);
    t_star[b] = (theta_star - theta_hat) / se_star;
  }
  
  // Sort t_star for shortest interval search
  NumericVector sorted_t = clone(t_star);
  std::sort(sorted_t.begin(), sorted_t.end());
  
  std::pair<double, double> si = shortest_interval(sorted_t, alpha);
  double lower = theta_hat - si.second * se_hat;
  double upper = theta_hat - si.first * se_hat;
  
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat short_bootstrap_t_grid(NumericVector x, double parm, int B = 100,
                           double alpha_start = 0.0005, double alpha_end = 0.20,
                           int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  double theta_hat = mean(x);
  double se_hat = sd(x) / sqrt(n);
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  
  // Bootstrap
  NumericVector t_star(B);
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      x_star[i] = x[idx];
    }
    double theta_star = mean(x_star);
    double se_star = sd(x_star) / sqrt(n);
    t_star[b] = (theta_star - theta_hat) / se_star;
  }
  
  // Sort t_star for shortest interval search
  NumericVector sorted_t = clone(t_star);
  std::sort(sorted_t.begin(), sorted_t.end());
  
  
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {

    // Shortest interval search
    std::pair<double, double> si = shortest_interval(sorted_t, alpha);
    double lower = theta_hat - si.second * se_hat;
    double upper = theta_hat - si.first * se_hat;
    
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// [[Rcpp::export]]
NumericVector short_percentile_bootstrap(NumericVector x, double parm, int B = 100, double alpha = 0.05) {
  int n = x.size();
  NumericVector theta_star(B);
  
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      x_star[i] = x[idx];
    }
    theta_star[b] = mean(x_star);
  }
  
  // Sort bootstrap estimates for shortest interval search
  NumericVector sorted_theta = clone(theta_star);
  std::sort(sorted_theta.begin(), sorted_theta.end());
  
  std::pair<double, double> si = shortest_interval(sorted_theta, alpha);
  double lower = si.first;
  double upper = si.second;
  
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat short_percentile_bootstrap_grid(NumericVector x, double parm, int B = 100,
                                 double alpha_start = 0.0005, double alpha_end = 0.20,
                                 int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data parameters
  int n = x.size();
  
  // Bootstrap
  NumericVector theta_star(B);
  
  for(int b = 0; b < B; ++b) {
    NumericVector x_star(n);
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));
      x_star[i] = x[idx];
    }
    theta_star[b] = mean(x_star);
  }
  
  // Sort bootstrap estimates for shortest interval search
  NumericVector sorted_theta = clone(theta_star);
  std::sort(sorted_theta.begin(), sorted_theta.end());
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    
    // Shortest interval search
    std::pair<double, double> si = shortest_interval(sorted_theta, alpha);
    double lower = si.first;
    double upper = si.second;
    
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// [[Rcpp::export]]
NumericVector w_cpp(NumericVector X, double gamma) {
  return pow((X + gamma) * sum(1 / (X + gamma)), -1);
}

// [[Rcpp::export]]
NumericVector w_prime_cpp(NumericVector X, double gamma) {
  NumericVector a_term = 1 / (X + gamma);
  NumericVector a_term_prime = 1 / pow(X + gamma, 2);
  double S = sum(1 / (X + gamma));
  double S_term_prime = sum(1 / pow((X + gamma), 2));
  
  return -a_term_prime / S + a_term * S_term_prime / pow(S, 2);
}

// [[Rcpp::export]]
bool valid_gamma_cpp(NumericVector X, double gamma) {
  LogicalVector shifted_pos = (X + gamma > 0);
  LogicalVector shifted_neg = (X + gamma < 0);
  
  return is_true(all(shifted_pos)) || is_true(all(shifted_neg));
}

// [[Rcpp::export]]
double f_fun_cpp(NumericVector X, double gamma, double r0) {
  if (!valid_gamma_cpp(X, gamma)) return NA_REAL;
  int n = X.size();
  
  return -2 * sum(log(n * w_cpp(X, gamma))) - r0;
}

// [[Rcpp::export]]
double f_prime_fun_cpp(NumericVector X, double gamma, double r0) {
  if (!valid_gamma_cpp(X, gamma)) return NA_REAL;
  
  return -2 * sum(w_prime_cpp(X, gamma) / w_cpp(X, gamma));
}


// [[Rcpp::export]]
DataFrame Newton_Fun_cpp(Function f, Function f_prime,
                         NumericVector gamma_start,
                         NumericVector X, double r0,
                         double epsilon = 1e-7, int max_iter = 100) {
  
  int m = gamma_start.size();
  NumericVector gamma_new_vec(m);
  
  for (int i = 0; i < m; i++) {
    double gamma_old = gamma_start[i];
    double gamma_new = NA_REAL;
    
    // Ensure initial guess is in valid region
    if (gamma_start[i] < -max(X)) {
      gamma_new = gamma_old - 100 * epsilon;
    } else if (gamma_start[i] > -min(X)) {
      gamma_new = gamma_old + 100 * epsilon;
    } else {
      Rcout << "Invalid choice of gamma_start; must be < -max(X) or > -min(X).\n";
      gamma_new_vec[i] = NA_REAL;
      continue;
    }
    
    int iter = 0;
    while (!NumericVector::is_na(gamma_new) &&
           !std::isinf(gamma_new) &&
           std::fabs(gamma_new - gamma_old) >= epsilon) {
      
      gamma_old = gamma_new;
      
      double fval = as<double>(f(X, gamma_old, r0));
      double fprime = as<double>(f_prime(X, gamma_old, r0));
      
      if (NumericVector::is_na(fval) || NumericVector::is_na(fprime)) {
        Rcout << "Invalid function value at iteration " << iter 
              << " for start index " << i << "\n";
        gamma_new = NA_REAL;
        break;
      }
      
      gamma_new = gamma_old - fval / fprime;
      
      iter++;
      if (iter > max_iter) {
        Rcout << "Max iterations reached for starting value " << i << "\n";
        gamma_new = NA_REAL;
        break;
      }
    }
    
    gamma_new_vec[i] = gamma_new;
  }
  
  return DataFrame::create(Named("Starting_Val") = gamma_start,
                           Named("Estimated_gamma") = gamma_new_vec);
}

// Define signatures
typedef double (*FunPtr)(NumericVector, double, double);

// Newton using C++ function pointers
DataFrame Newton_Fun_cpp_internal(FunPtr f, FunPtr f_prime,
                                  NumericVector gamma_start,
                                  NumericVector X, double r0,
                                  double epsilon = 1e-7, int max_iter = 100) {
  
  int m = gamma_start.size();
  NumericVector gamma_new_vec(m);
  
  for (int i = 0; i < m; i++) {
    double gamma_old = gamma_start[i];
    double gamma_new = NA_REAL;
    
    if (gamma_start[i] < -max(X)) {
      gamma_new = gamma_old - 100 * epsilon;
    } else if (gamma_start[i] > -min(X)) {
      gamma_new = gamma_old + 100 * epsilon;
    } else {
      gamma_new_vec[i] = NA_REAL;
      continue;
    }
    
    int iter = 0;
    while (!NumericVector::is_na(gamma_new) &&
           !std::isinf(gamma_new) &&
           std::fabs(gamma_new - gamma_old) >= epsilon) {
      
      gamma_old = gamma_new;
      
      double fval = f(X, gamma_old, r0);
      double fprime = f_prime(X, gamma_old, r0);
      
      if (NumericVector::is_na(fval) || NumericVector::is_na(fprime)) {
        gamma_new = NA_REAL;
        break;
      }
      
      gamma_new = gamma_old - fval / fprime;
      
      iter++;
      if (iter > max_iter) {
        gamma_new = NA_REAL;
        break;
      }
    }
    
    gamma_new_vec[i] = gamma_new;
  }
  
  return DataFrame::create(Named("Starting_Val") = gamma_start,
                           Named("Estimated_gamma") = gamma_new_vec);
}


// [[Rcpp::export]]
NumericVector emp_lik(NumericVector x, double parm, std::string type, double alpha = 0.05) {
  int n = x.size();
  double r0;
  if (type == "C" || type == "c" || type == "Chi" || type == "chi") {
    // Calling the C-level API (from Rmath.h), and those functions have slightly different signatures:
    // double qnchisq(double p, double df, double ncp, int lower_tail, int log_p);
    // double qnf(double p, double df1, double df2, double ncp, int lower_tail, int log_p);
    r0 = R::qchisq(1 - alpha, 1, 1, 0);
    // Rcout << r0 << "\n";
  } else if (type == "F" || type == "f") {
    r0 = R::qf(1 - alpha, 1, n - 1, 1, 0);
    // Rcout << r0 << "\n";
  } else {
    throw std::runtime_error("Empirical likelihood critical value not supported");
  }
  
  double x_min = min(x);
  double x_max = max(x);
  double x_range = x_max - x_min;
  
  NumericVector gamma_starts(2);
  
  gamma_starts[0] = -x_max - 0.0001 * x_range;
  gamma_starts[1] = -x_min + 0.0001 * x_range;
  
  
  DataFrame gammaVals = Newton_Fun_cpp_internal(f_fun_cpp, f_prime_fun_cpp, gamma_starts, x, r0);
 
  NumericVector est_gamma = gammaVals["Estimated_gamma"];
  double lower = sum(w_cpp(x, est_gamma[1]) * x);
  double upper = sum(w_cpp(x, est_gamma[0]) * x);
  double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
  double length = upper - lower;
  
  return NumericVector::create(lower, upper, captured, length);
}

// [[Rcpp::export]]
arma::mat emp_lik_grid(NumericVector x, double parm, std::string type,
                       double alpha_start = 0.0005, double alpha_end = 0.20,
                       int alpha_length = 400) {
  // Output // matrix: rows = coefficients, cols = output fields
  arma::mat alpha_grid(alpha_length, 5);
  
  // Data stats
  int n = x.size();
  double x_min = min(x);
  double x_max = max(x);
  double x_range = x_max - x_min;
  
  NumericVector gamma_starts(2);
  
  gamma_starts[0] = -x_max - 0.0001 * x_range;
  gamma_starts[1] = -x_min + 0.0001 * x_range;
  
  // Initialize critical value
  double r0;
  
  // Step size
  double alpha_step = (alpha_end - alpha_start) / (alpha_length - 1);
  
  // Integer index for rows
  int row = 0; 
  for (double alpha = alpha_start; alpha <= alpha_end + 1e-12; alpha += alpha_step) {
    // Which EL are we doing?
    if (type == "C" || type == "c" || type == "Chi" || type == "chi") {
      r0 = R::qchisq(1 - alpha, 1, 1, 0);
    } else if (type == "F" || type == "f") {
      r0 = R::qf(1 - alpha, 1, n - 1, 1, 0);
    } else {
      throw std::runtime_error("Empirical likelihood critical value not supported");
    }
    
    // Return Newton
    DataFrame gammaVals = Newton_Fun_cpp_internal(f_fun_cpp, f_prime_fun_cpp, gamma_starts, x, r0);
    
    NumericVector est_gamma = gammaVals["Estimated_gamma"];
    
    double lower = sum(w_cpp(x, est_gamma[1]) * x);
    double upper = sum(w_cpp(x, est_gamma[0]) * x);
    double captured = (lower < parm && parm < upper) ? 1.0 : 0.0;
    double length = upper - lower;
    
    // Store in output matrix
    alpha_grid(row, 0) = lower;
    alpha_grid(row, 1) = upper;
    alpha_grid(row, 2) = captured;
    alpha_grid(row, 3) = length;
    alpha_grid(row, 4) = alpha;
    
    row++;
  }
  
  return alpha_grid;
}

// Testing
/*** R
# Simulations to explore CI performance
####

## Data-generating parameters
n <- 3
d <- 1

## Hypothesis Testing
nsim <- 1000
alpha <- 0.05
mu0_vec <- 1:10
r0 <- qchisq(1 - alpha, df = d)

## True parameter value
mu0 <- 4

samp <- rgamma(n, shape = 2, scale = 2)

emp_lik(samp, mu0, type = "F")
*/
