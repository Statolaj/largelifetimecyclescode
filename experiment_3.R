####### Initialization ####### 

library(rstudioapi)
library(TDA)
library(spatstat)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_file.R")


####### Simulated correlation #######

n_values <- seq(100, 1000, by = 100)
average_correlations <- numeric(length(n_values))
m <- 3

for (j in seq_along(n_values)) {
  print(j)
  n <- n_values[j]
  K <- 20
  correlations <- numeric(K)
  
  # Bound
  beta_low <- 0.5 + 1/(2*m)
  bound <- n^(-beta_low)
  bound^2
  for (k in 1:K) {
    poisson_samples <- replicate(100, rpoispp(10 * n, win = square(1)))
    
    # Compute the correlation between the death time of the largest feature and the death time of the second largest feature
    largest_deathtimes <- numeric(100)
    second_largest_deathtimes <- numeric(100)
    i=1
    for (i in 1:100) {
      points <- cbind(poisson_samples[,i]$x, poisson_samples[,i]$y)
      result <- maximal_lifetime_2(points, bound, m)
      if (length(result$largest_lifetimes) == 2) {
        largest_deathtimes[i] <- result$corresponding_deathtimes[1]
        second_largest_deathtimes[i] <- result$corresponding_deathtimes[2]
      } else {
        largest_deathtimes[i] <- NA
        second_largest_deathtimes[i] <- NA
      }
    }
    valid_indices <- which(!is.na(largest_deathtimes) & !is.na(second_largest_deathtimes))
    if (length(valid_indices) > 1) {
      correlations[k] <- cor(largest_deathtimes[valid_indices], second_largest_deathtimes[valid_indices])
    } else {
      correlations[k] <- NA
    }
  }
  
  # Calculate average correlation for the current n
  average_correlations[j] <- mean(correlations, na.rm = TRUE)
}

# Plot the average correlation against the value of n
plot(n_values, average_correlations, type = "b", xlab = "Value of n", ylab = "Average Correlation", 
     main = " Correlation of largest and second largest deathtime", ylim = c(-0.1, max(average_correlations)))
