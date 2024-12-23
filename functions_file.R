####### Conjectured upper bound for maximal lifetime ####### 

conjectured_lmax <- function(m) {
  # Generate m equally spaced angles in the range [0, 2π]
  angles <- seq(0, 2 * pi, length.out = m + 1)[1:m]
  # Calculate the x and y coordinates using sine and cosine
  x <- cos(angles)
  y <- sin(angles)
  # Compute the Euclidean distance between the first two points
  dist <- dist(rbind(c(x[1], y[1]), c(x[2], y[2])))
  # Return half of the distance
  return(1 - dist / 2)
}


###### Generate points uniformly inside a unit circle #######

generate_points_in_circle <- function(num_points = 4) {
  angles <- runif(num_points, 0, 2 * pi)
  radii <- sqrt(runif(num_points))
  x <- radii * cos(angles)
  y <- radii * sin(angles)
  return(cbind(x, y))
}

scale_points_to_max_distance <- function(flattened_points) {
  points <- matrix(flattened_points, ncol = 2, byrow = TRUE)
  distances <- sqrt(rowSums(points^2))
  max_distance <- max(distances)
  # If the max distance is 0, return the points unchanged (all points at origin)
  if (max_distance == 0) {
    return(flattened_points)
  }
  scaled_points <- points / max_distance
  return(as.vector(t(scaled_points)))
}


####### Plot points inside circle ####### 

plot_points_in_circle <- function(points, title) {
  # Create the plot
  plot(NA, xlim = c(-1, 1), ylim = c(-1.1, 1.1), asp = 1,
       xlab = "X-axis", ylab = "Y-axis",
       main = title)
  
  # Draw the unit circle
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, lty = 2, fg = "black")
  
  # Plot the points
  points(points[,1], points[,2], col = "red", pch = 19)
  
  # Add grid lines
  grid()
}


####### Extract largest lifetimes for unit circle ####### 

simple_maximal_lifetime_cech <- function(flattened_points) {
  # Reshape the flattened vector back into a matrix of points
  points <- matrix(flattened_points, ncol = 2, byrow = TRUE)
  
  # Calculate the distance of each point from the origin
  distances <- sqrt(rowSums(points^2))
  
  # Define a large penalty
  penalty <- 1e6
  
  # Check if any points are outside the unit circle
  if (any(distances > 1)) {
    return(-penalty)
  }
  
  # Proceed with the original logic if all points are within the unit circle
  diag <- alphaComplexDiag(points, library = "GUDHI", printProgress = FALSE)
  #diag <- ripsDiag(points, library = "GUDHI", maxdimension = 1, maxscale = 1, printProgress = FALSE)
  holes <- diag[["diagram"]][diag[["diagram"]][, 1] == 1, ]
  holes <- matrix(data = holes, ncol = 3)
  
  if (nrow(holes) > 0) {
    lifetimes <- sqrt(holes[, 3]) - sqrt(holes[, 2])
    return(max(lifetimes))
  } else {
    return(0)
  }
}


####### Extract two largest lifetimes ####### 

maximal_lifetime_2 <- function(flattened_points, bound, m) {
  # Reshape the flattened vector back into a matrix of points
  diag <- alphaComplexDiag(matrix(flattened_points, ncol = 2, byrow = TRUE),
                           library = c("GUDHI", "Dionysus"), location = TRUE,
                           printProgress = FALSE)
  
  # Extract holes of dimension 1 and keep their locations
  holes <- diag[["diagram"]][diag[["diagram"]][, 1] == 1, ]
  locations <- diag[["cycleLocation"]][diag[["diagram"]][, 1] == 1]
  
  # Check if holes is NULL or has zero rows
  if (is.null(holes) || length(holes) == 0 || nrow(holes) == 0) {
    return(list(largest_deathtimes = numeric(0), corresponding_deathtimes = numeric(0)))
  }
  
  # Remove holes consisting of more than m points
  holes <- holes[sapply(locations, nrow) <= m, ]
  
  # Filter out holes with death time greater than the specified bound
  holes <- holes[sqrt(holes[, 3]) < bound, ]
  
  # Check if there are any remaining holes
  if (nrow(holes) > 0) {
    # Calculate the lifetimes
    lifetimes <- sqrt(holes[, 3]) - sqrt(holes[, 2])
    
    # Get the two largest death times
    if (length(lifetimes) >= 2) {
      indices <- order(lifetimes, decreasing = TRUE)[1:2]
    } else {
      indices <- order(lifetimes, decreasing = TRUE)[1]
    }
    
    largest_lifetimes <- lifetimes[indices]
    corresponding_deathtimes <- sqrt(holes[indices, 3])
    
    return(list(largest_lifetimes = largest_lifetimes, corresponding_deathtimes = corresponding_deathtimes))
  } else {
    return(list(largest_lifetimes = numeric(0), corresponding_deathtimes = numeric(0)))
  }
}