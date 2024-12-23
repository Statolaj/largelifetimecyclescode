####### Initialization ####### 

library(rstudioapi)
library(TDA)
library(spatstat)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_file.R")

####### Part 1 #######

# Define the intensities, t_vals and number of repetitions
n_values <- c(20, 50, 100, 200, 500, 1000)
num_reps <- 5000
t_values <- seq(0.02, 0.1, by = 0.001)  

# For loop
par(mfrow = c(2, 3))
for (n in n_values) {
  cat("Processing n =", n, "\n")
  
  # Threshold for filtering holes
  death_threshold <- n^(-0.75)
  lmax = 1-sqrt(3)/2
  
  replicate_stats <- list() # store results per replicate
  
  for (rep_i in 1:num_reps) {
    # Generate a Poisson point process with intensity n on [0,1]^2
    pp <- rpoispp(lambda = n, win = owin(c(0,1), c(0,1)))
    coords <- cbind(pp$x, pp$y)
    
    # Compute alpha complex persistence diagram up to dimension 1
    # (If a true Čech filtration is required and available, substitute here)
    diag <- alphaComplexDiag(X = coords, maxdimension = 1)
    dgm <- diag$diagram
    
    # Extract dimension 1 features (holes)
    dgm_1 <- dgm[dgm[,1] == 1, , drop = FALSE]
    dgm_1
    
    # Filter holes: keep only those with death time < n^{-0.67}
    filtered_holes <- dgm_1[sqrt(dgm_1[,3]) < death_threshold, , drop = FALSE]
    
    # Here, we need to define what corresponds to the "largest additive lifetime".
    # Suppose "lifetime" is death-birth. Let's assume the largest lifetime is max(death-birth) among filtered holes.
    # You may need to adjust this according to your interpretation of L^+_{max} - L^(1)_n.
    if (nrow(filtered_holes) > 0) {
      lifetimes <- (sqrt(filtered_holes[,3]) - sqrt(filtered_holes[,2]))/death_threshold
      largest_lifetime <- max(lifetimes)
    } else {
      largest_lifetime <- 0
    }
    
    # Store the largest lifetime of this replicate
    if (largest_lifetime < lmax){
      replicate_stats[[rep_i]] <- largest_lifetime
    }
  }
  
  # Once we have 30 samples of the largest lifetime for given n:
  largest_lifetimes <- unlist(replicate_stats)
  
  # Estimate S(t) = P(L^+_{max} - L^(1)_n > t)
  # For each t, we compute the empirical survival probability
  survival_probs <- sapply(t_values, function(t) mean((lmax - largest_lifetimes) > t))
  
  # Compute x and y
  x <- log(t_values)
  y <- log(-log(survival_probs))
  
  # Filter out infinite y-values
  valid_indices <- is.finite(y)  # Check for finite values
  x_filtered <- x[valid_indices]
  y_filtered <- y[valid_indices]
  
  fit <- lm(y_filtered ~ x_filtered)
  slope <- round(coef(fit)[2],3)
  
  # Plot the filtered data
  plot(x_filtered, y_filtered, xlim = c(-4,-2), ylim=c(-4,0), main = paste("Intensity n =", n,"and slope q =", slope), 
       xlab = "log(t)", ylab = "log(-log(S(t)))", cex = 0.75, cex.lab = 1.15)
  
  # Add a linear regression line
  fit <- lm(y_filtered ~ x_filtered)
  abline(fit, col = "blue", lty = 1, lwd = 1.3)
}

####### Part 2 #######


# Define the intensities and number of repetitions
n_values <- c(1000,2000,3500,5000)
num_reps <- 10000

# For loop
par(mfrow = c(3, 4))
for (m in c(4,5,6)) {
for (n in n_values) {
  cat("Processing n =", n, "\n")
  
  # Threshold for filtering holes
  if (m == 4){
  t_values <- seq(0.05, 0.17, by = 0.002)  
  death_threshold <- n^(-0.66)
  lmax = 0.293
  }
  if (m == 5){
    t_values <- seq(0.10, 0.22, by = 0.002)  
    death_threshold <- n^(-0.62)
    lmax = 0.412
  }
  if (m == 6){
    t_values <- seq(0.15, 0.3, by = 0.002)  
    death_threshold <- n^(-0.6)
    lmax = 0.500
  }
  
  replicate_stats <- list() # store results per replicate
  
  for (rep_i in 1:num_reps) {
    # Generate a Poisson point process with intensity n on [0,1]^2
    pp <- rpoispp(lambda = n, win = owin(c(0,1), c(0,1)))
    coords <- cbind(pp$x, pp$y)
    
    # Compute alpha complex persistence diagram up to dimension 1
    # (If a true Čech filtration is required and available, substitute here)
    diag <- alphaComplexDiag(X = coords, maxdimension = 1)
    dgm <- diag$diagram
    
    # Extract dimension 1 features (holes)
    dgm_1 <- dgm[dgm[,1] == 1, , drop = FALSE]
    dgm_1
    
    # Filter holes: keep only those with death time < n^{-0.67}
    filtered_holes <- dgm_1[sqrt(dgm_1[,3]) < death_threshold, , drop = FALSE]
    
    # Here, we need to define what corresponds to the "largest additive lifetime".
    # Suppose "lifetime" is death-birth. Let's assume the largest lifetime is max(death-birth) among filtered holes.
    # You may need to adjust this according to your interpretation of L^+_{max} - L^(1)_n.
    if (nrow(filtered_holes) > 0) {
      lifetimes <- (sqrt(filtered_holes[,3]) - sqrt(filtered_holes[,2]))/death_threshold
      largest_lifetime <- max(lifetimes)
    } else {
      largest_lifetime <- 0
    }
    
    # Store the largest lifetime of this replicate
    if (largest_lifetime < lmax){
      replicate_stats[[rep_i]] <- largest_lifetime
    }
  }
  
  # Once we have 30 samples of the largest lifetime for given n:
  largest_lifetimes <- unlist(replicate_stats)
  
  # Estimate S(t) = P(L^+_{max} - L^(1)_n > t)
  # For each t, we compute the empirical survival probability
  survival_probs <- sapply(t_values, function(t) mean((lmax - largest_lifetimes) > t))
  
  # Compute x and y
  x <- log(t_values)
  y <- log(-log(survival_probs))
  
  # Filter out infinite y-values
  valid_indices <- is.finite(y)  # Check for finite values
  x_filtered <- x[valid_indices]
  y_filtered <- y[valid_indices]
  
  fit <- lm(y_filtered ~ x_filtered)
  slope <- round(coef(fit)[2],2)
  
  # Plot the filtered data
  plot(x_filtered, y_filtered, xlim = c(-3.5,-1), ylim=c(-6.5,0), main = paste("m =", m, ", n =", n,", q =", slope), 
         xlab = "log(t)", ylab = "log(-log(S(t)))", cex = 0.75, cex.lab = 1.2)

  fit <- lm(y_filtered ~ x_filtered)
  abline(fit, col = "blue", lty = 1, lwd = 1.3)
 }
}
