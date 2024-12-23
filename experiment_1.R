####### Initialization ####### 

library(rstudioapi)
library(TDA)
library(optimization)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_file.R")


####### Simulated annealing #######

m <- 6
points <- generate_points_in_circle(m)
initial_guess <- as.vector(points)
lower_bounds <- rep(-1, length(initial_guess))
upper_bounds <- rep(1, length(initial_guess))

sim <- optim_sa(simple_maximal_lifetime_cech, start = initial_guess, maximization = TRUE, trace = FALSE,lower = lower_bounds,upper = upper_bounds,
                control = list(ac_acc = 0.00001, t0 = 500, r=0.95, nlimit = 30000, maxgood = 4000))


distances <- sqrt(rowSums(matrix(sim$par, ncol = 2, byrow = TRUE)^2))
max_distance <- max(distances)
print(sim$function_value/max_distance)


####### Generating plots #######

#par(mfrow = c(1, 4))
plot_points_in_circle(points, paste("m =", m))
points(matrix(scale_points_to_max_distance(as.vector(sim$par)), ncol = 2, byrow = T), col = "blue", pch = 16)
