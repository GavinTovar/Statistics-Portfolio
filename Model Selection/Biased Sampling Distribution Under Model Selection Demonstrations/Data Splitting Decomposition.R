# Decomposing the 4th exploration
library(mvtnorm)
library(ggplot2)
library(reshape2)
# Parameter settings
set.seed(10262025)
mean.vec <- c(1, 1, 1)
Sigma.Mat <- matrix(0.5, nrow = 3, ncol = 3)
diag(Sigma.Mat) <- 1
n <- 50
epsilon.var <- 25
Beta.seq <- c(1, -3, -1, -2)

mu.vec = mean.vec 
cov.mat = Sigma.Mat
epsilon.error = epsilon.var 
beta.seq = Beta.seq
MS.type = "forward"
nsim = 10000
criteria = "AIC"

# Split percentage
DSP <- 0.50

### Empty object creation
# Object for chosen covariates representing the model selected
Reg.Sel.List <- rep(NA, nsim)

# Object for beta coefficients and t-stats for all obervations and models nsim x 4 x 8
Coef.arr.train <- Tval.arr.train <- array(NA, dim = c(nsim, length(beta.seq), 8))
Coef.arr.test <- Tval.arr.test <- array(NA, dim = c(nsim, length(beta.seq), 8))

### Start Simulations
for(i in 1:nsim){
  ## Random MVN Sample
  x.mat <- rmvnorm(n = n, mean = mu.vec, sigma = cov.mat)
  
  # ABC labeling of regressors unless there is a lot
  if(length(mu.vec) > 24){
    colnames(x.mat) <- paste("X", 1:length(mu.vec), sep = "")
  }
  colnames(x.mat) <- LETTERS[1:length(mu.vec)]
  
  # Creating response
  y.vec <- beta.seq[1] + beta.seq[-1]%*%t(x.mat) + rnorm(n, mean = 0, sd = sqrt(epsilon.error))
  
  # Creating dataframe
  da.data <- cbind.data.frame(Y = as.vector(y.vec), x.mat)
  
  # Split data into training and test sets
  index <- sample(1:n, size = n * DSP)
  train.data <- da.data[index, ]
  test.data <- da.data[-index, ]
  
  # Performing specified model selection technique
  if(MS.type %in% c("forward", "backward", "both")){
    # Create both null models
    train.null.model <- lm(Y ~ 1, data = train.data)
    test.null.model <- lm(Y ~ 1, data = test.data)
    
    # Only use training data/model for model selection
    Mod.Select <- step(object = train.null.model, scope = as.formula(paste("~", paste(colnames(x.mat), collapse = " + "))), 
                       direction = MS.type, trace = 0, k = ifelse(criteria == "AIC", 2, log(n)))
    
    # Noting which regressors were chosen in the selected mode
    Reg.Sel.List[i] <- paste(sort(attr(Mod.Select$terms, "term.labels")), collapse = "")
    
    
    # List of possible models
    model.list <- list(c(), c("A"), c("B"), c("C"), 
                       c("A", "B"), c("A", "C"), c("B", "C"), 
                       c("A", "B", "C"))
    
    # Keeping coefficents from null models
    # Train
    Coef.arr.train[i, 1, 1] <- summary(train.null.model)$coefficients[,1]
    Tval.arr.train[i, 1, 1] <- summary(train.null.model)$coefficients[,3]
    # Test
    Coef.arr.test[i, 1, 1] <- summary(test.null.model)$coefficients[,1]
    Tval.arr.test[i, 1, 1] <- summary(test.null.model)$coefficients[,3]
    
    for(j in 2:8){
      Cur.Regs <- model.list[[j]]
      # Train
      train.lm.curr <- summary(lm(formula = as.formula(paste("Y ~", paste(Cur.Regs, collapse = " + "))), data = train.data))
      train.indd <- match(attr(train.lm.curr$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics
      Coef.arr.train[i, c(1, train.indd+1), j] <- train.lm.curr$coefficients[,1]
      Tval.arr.train[i, c(1, train.indd+1), j] <- train.lm.curr$coefficients[,3]
      
      # Test
      test.lm.curr <- summary(lm(formula = as.formula(paste("Y ~", paste(Cur.Regs, collapse = " + "))), data = test.data))
      test.indd <- match(attr(test.lm.curr$terms, "term.labels"), colnames(x.mat))
      
      # Saving relevant statistics
      Coef.arr.test[i, c(1, test.indd+1), j] <- test.lm.curr$coefficients[,1]
      Tval.arr.test[i, c(1, test.indd+1), j] <- test.lm.curr$coefficients[,3]
    }
  }
  # Print statement
  if(i %% 500 == 0) print(i)
}

# Labeling coefficient array
dimnames(Coef.arr.test) <- dimnames(Coef.arr.train) <- list(1:nsim, 
                                                            c("Int", "A", "B", "C"), 
                                                            c("1", "A", "B", "C", "AB", "AC", "BC", "ABC"))
# Average estimated coefficients without model selection
# Train
mat1 <- t(round(apply(Coef.arr.train, c(2,3), mean), 2))

# Test
mat2 <- t(round(apply(Coef.arr.test, c(2,3), mean), 2))

# How often was each model chosen?
table(Reg.Sel.List)

# Possible Models & Coeffs
Mods <- c("1", "A", "B", "C", "AB", "AC", "BC", "ABC")
Coeffs <- c("Int", "A", "B", "C")

# Noting which models were chosen via model selection
ind.1 <- which(Reg.Sel.List == "")
ind.A <- which(Reg.Sel.List == "A")
ind.B <- which(Reg.Sel.List == "B")
ind.C <- which(Reg.Sel.List == "C")
ind.AB <- which(Reg.Sel.List == "AB")
ind.AC <- which(Reg.Sel.List == "AC")
ind.BC <- which(Reg.Sel.List == "BC")
ind.ABC <- which(Reg.Sel.List == "ABC")

# Averaging coefficients for models chosen post-model selection
# Train
mat3 <- rbind(
  t(round(apply(Coef.arr.train[ind.1,,1], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.A,,2], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.B,,3], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.C,,4], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.AB,,5], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.AC,,6], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.BC,,7], 2, mean), 2)),
  t(round(apply(Coef.arr.train[ind.ABC,,8], 2, mean), 2)))
rownames(mat3) <- Mods
colnames(mat3) <- Coeffs

# Test
mat4 <- rbind(
  t(round(apply(Coef.arr.test[ind.1,,1], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.A,,2], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.B,,3], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.C,,4], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.AB,,5], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.AC,,6], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.BC,,7], 2, mean), 2)),
  t(round(apply(Coef.arr.test[ind.ABC,,8], 2, mean), 2)))
rownames(mat4) <- Mods
colnames(mat4) <- Coeffs

# Combine matrices
gbg <- list(Train.Pre = mat1, Test.Pre = mat2,
            Train.Post = mat3, Test.Post = mat4)
gbg

Plot.Function <- function(array1, array2, col = NULL){
  # Helper function to process list of vectors from indexed arrays
  process.vectors <- function(array){
    ### Creating post-MS coef matricies
    vectors <- list(array[ind.1, col, 1],
                    array[ind.A, col, 2],
                    array[ind.B, col, 3],
                    array[ind.C, col, 4],
                    array[ind.AB, col, 5],
                    array[ind.AC, col, 6],
                    array[ind.BC, col, 7],
                    array[ind.ABC, col, 8])
    
    # Determine the maximum length of the vectors
    max_len <- max(sapply(vectors, length))
    
    # Pad shorter vectors with NA to match the maximum length
    padded_vectors <- lapply(vectors, function(x) {
      length(x) <- max_len
      return(x)
    })
    
    # Combine the padded vectors into a matrix with 8 columns
    matrix_data <- do.call(cbind, padded_vectors)
    
    #
    colnames(matrix_data) <- Mods
    
    # Return resulting matrix
    return(matrix_data)
  }
  
  matrix1 <- process.vectors(array1)
  matrix2 <- process.vectors(array2)
  
  # Helper function to process each matrix
  process_matrix <- function(matrix, dataset_name) {
    CoeffMat <- matrix
    
    # Remove NA columns
    CoeffMat <- CoeffMat[, colSums(is.na(CoeffMat)) != nrow(CoeffMat)]
    
    # Reshape data into long format
    data_long <- reshape2::melt(CoeffMat, variable.name = "Facet", value.name = "Value")[,-1]
    
    # Turn col input into correct name
    CoefName <- switch(col, "Intercept", "Beta1", "Beta2", "Beta3")
    
    colnames(data_long) <- c("Model", CoefName)
    
    # Add a column to indicate which dataset the data belongs to
    data_long$Dataset <- dataset_name
    
    # Remove rows w/ padded NAs
    data_long <- na.omit(data_long)
    
    return(data_long)
  }
  
  # Helper function to process each array
  process_array <- function(array, dataset_name) {
    CoeffMat <- array[, col, ]
    
    # Remove NA columns
    CoeffMat <- CoeffMat[, colSums(is.na(CoeffMat)) != nrow(CoeffMat)]
    
    # Reshape data into long format
    data_long <- reshape2::melt(CoeffMat, variable.name = "Facet", value.name = "Value")[,-1]
    
    # Turn col input into correct name
    CoefName <- switch(col, "Intercept", "Beta1", "Beta2", "Beta3")
    
    colnames(data_long) <- c("Model", CoefName)
    
    # Add a column to indicate which dataset the data belongs to
    data_long$Dataset <- dataset_name
    
    return(data_long)
  }
  
  # Process both arrays
  data1 <- process_array(array1, "Training Pre-MS")
  data2 <- process_array(array2, "Test Pre-MS")
  data3 <- process_matrix(matrix1, "Training Post-MS")
  data4 <- process_matrix(matrix2, "Test Post-MS")
  
  # Combine the datasets
  combined_data <- rbind(data1, data2, data3, data4)
  
  # Establish CoeffName in higher scope
  CoefName <- colnames(combined_data)[2]
  
  # Check the combined data structure
  print(head(combined_data))
  print(tail(combined_data))
        
  # Create the plot
  plot <- ggplot(combined_data, aes(x = !!sym(CoefName), fill = Dataset, color = Dataset)) +
    geom_density(alpha = 0.1333) +
    facet_wrap(~ Model, labeller = labeller(Model = function(x) paste("Model", x))) +
    labs(title = "Facets of Each Different Model, Overlaid Densities for Each Dataset",
         x = paste(CoefName, "Coefficients"),
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(plot)
}

# Intercept
Plot.Function(Coef.arr.train, Coef.arr.test, col = 1)
# Beta1
Plot.Function(Coef.arr.train, Coef.arr.test, col = 2)
# Beta2
Plot.Function(Coef.arr.train, Coef.arr.test, col = 3)
# Beta3
Plot.Function(Coef.arr.train, Coef.arr.test, col = 4)
