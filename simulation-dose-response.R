rm(list = ls())
set.seed(123)
library(tidyverse)

# simulation ----------------
X1 <- c(19.3, 19.3, 16.3, 13.2, 22.7, 15.3, 11.9, 16.9, 17.2, 31.4, 19.6, 15.8,
        13.5, 18.4, 17.9, 22.6, 19.5, 17.3, 20.9, 28.0, 27.5, 23.9, 19.0, 17.1)
X2 <- c(210.0, 91.2, 372.8, 228.5, 112.5, 175.3, 185.8, 241.7, 205.1, 268.2, 120.9, 113.4,
        158.0, 196.4, 182.9, 333.3, 390.9, 113.3, 75.5, 100.7, 385.0, 416.3, 148.5, 153.0)
population <- c(12747, 10179, 17702, 6471, 12936, 22085, 19587, 16026,
                13132, 8127, 14983, 18947, 28541, 15101, 8976, 10474,
                4031, 21185, 9726, 10780, 8348, 9971, 6677, 4950)/10000
Y0 <- round(X1 * population)
Z0 <- round(X2 * population)
sim_size <- 10000
gamma <- log(0.6)
eta <- ksi <- 0.1 # rates of imprecise tests
alpha <- rbeta(n = 24, shape1 = 0.5, shape2 = 0.5)

simulated_data_DR <- map(1:sim_size, function(j){
  D0 <- runif(24, 0.01, 0.4)
  D1 <- runif(24, 0.6, 1)
  sim_Y00 <- rmultinom(1, size = sum(Y0), prob = Y0/sum(Y0)) * (2 * population) 
  sim_Yd0 <- sim_Y00 * exp(gamma * D0)
  sim_Zd0 <- rmultinom(1, size = sum(Z0), prob = Z0/sum(Z0)) / (2 * population)
  sim_Zd1 <- alpha * sim_Zd0
  sim_Yd1 <- alpha * sim_Y00 * exp(gamma * D1)
  obs_Yd0 <- (1-ksi) * sim_Yd0 + eta * sim_Zd0
  obs_Yd1 <- (1-ksi) * sim_Yd1 + eta * sim_Zd1
  obs_Zd0 <- ksi * sim_Yd0 + (1-eta) * sim_Zd0
  obs_Zd1 <- ksi * sim_Yd1 + (1-eta) * sim_Zd1
  A <- sample(rep(0:1, each = 12), size = 24, replace = F)
  data.frame(A = A,
             D = A * D1 + (1-A) * D0,
             Y = A * obs_Yd1 + (1-A) * obs_Yd0,
             Z = A * obs_Zd1 + (1-A) * obs_Zd0)
})

dose_response <- function(Y, Z, D, A, X = NULL){
  if(min(Y) < 1e-5) {return(rep(NA, 5))}
  gamma_grid <- seq(-5, 5, by = 0.02)
  L <- log(Y/Z)
  m <- length(A)
  m1 <- sum(A)
  m0 <- m - m1
  if(is.null(X)){
    p_value <- map_dbl(gamma_grid, function(gamma){
      est <- mean(L[A==1] - gamma * D[A==1]) - mean(L[A==0] - gamma * D[A==0])
      var <- var(L[A==1] - gamma * D[A==1])/m1 + var(L[A==0] - gamma * D[A==0])/m0
      2 * pt(-abs(est)/sqrt(var), df = m)
    })
  } else {
    p_value <- map_dbl(gamma_grid, function(gamma){
      X <- as.matrix(X)
      V_X <- cov(X)
      beta <- solve(V_X) %*% cov(X, L - gamma * D)
      est <- mean(L[A==1] - gamma * D[A==1] - X[A==1,] %*% beta) - 
        mean(L[A==0] - gamma * D[A==0] - X[A==0,] %*% beta)
      var <- var(L[A==1] - gamma * D[A==1]- X[A==1,] %*% beta)/m1 * (m1-1)/(m1-1-ncol(X)) + 
        var(L[A==0] - gamma * D[A==0] - X[A==0,] %*% beta)/m0 * (m0-1)/(m0-1-ncol(X))
      2 * pt(-abs(est)/sqrt(var), df = m)
    })
  }
  est <- gamma_grid[which.max(p_value)]
  CI.lower <- range(gamma_grid[p_value >= 0.05])[1]
  CI.upper <- range(gamma_grid[p_value >= 0.05])[2]
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= gamma & CI.upper >= gamma
  c(est - gamma, CI.lower, CI.upper, power, coverage)
}

result <- map_dfc(simulated_data_DR, ~dose_response(.$Y, .$Z, .$D, .$A, cbind(population)))
map_dbl((1:nrow(result)), ~mean(unlist(result[.,]), na.rm = T))
result <- map_dfc(simulated_data_DR, ~dose_response(.$Y, .$Z, .$D, .$A))
map_dbl((1:nrow(result)), ~mean(unlist(result[.,]), na.rm = T))

# Y <- simulated_data_DR[[1]]$Y
# Z <- simulated_data_DR[[1]]$Z
# A <- simulated_data_DR[[1]]$A
# D <- simulated_data_DR[[1]]$D
# X <- cbind(population)
# dose_response(Y,Z,D,A,X)

