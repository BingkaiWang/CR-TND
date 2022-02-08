rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack)
library(lme4)
d <- readxl::read_xlsx("AWED data by cluster_BingkaiWang.xlsx", sheet = 3)
d <- d %>% group_by(cluster) %>% 
  summarise(intervention = mean(intervention),
            wei = mean(wei_activity_cat, na.rm=T)/5 + 0.1,
            test_neg = sum(test_neg),
            test_pos = sum(test_pos)) %>%
  as.data.frame()
d[1:24,]
X1 <- c(19.3, 19.3, 16.3, 13.2, 22.7, 15.3, 11.9, 16.9, 17.2, 31.4, 19.6, 15.8,
        13.5, 18.4, 17.9, 22.6, 19.5, 17.3, 20.9, 28.0, 27.5, 23.9, 19.0, 17.1)
X2 <- c(210.0, 91.2, 372.8, 228.5, 112.5, 175.3, 185.8, 241.7, 205.1, 268.2, 120.9, 113.4,
        158.0, 196.4, 182.9, 333.3, 390.9, 113.3, 75.5, 100.7, 385.0, 416.3, 148.5, 153.0)
population <- c(12747, 10179, 17702, 6471, 12936, 22085, 19587, 16026,
                13132, 8127, 14983, 18947, 28541, 15101, 8976, 10474,
                4031, 21185, 9726, 10780, 8348, 9971, 6677, 4950)
prop_children <- c(21, 21, 21, 20, 21, 20, 21, 21, 21, 22, 21, 20, 
                   21, 21, 22, 22, 23, 21, 22, 22, 23, 23, 24, 24)

# estimators ----------------
odds_ratio_est <- function(Y, Z, A, lambda = 1){
  m <- length(Y)
  m1 <- sum(A)
  m0 <- m - m1
  Y1 <- Y[A==1]
  Y0 <- Y[A==0]
  Z1 <- Z[A==1]
  Z0 <- Z[A==0]
  est <- log((sum(Y1)/sum(Y0))/(sum(Z1)/sum(Z0)))
  # calculate variance below
  n_Y <- sum(Y)
  n_Z <- sum(Z)
  V_Y <- var(Y1)*m1/m + var(Y0)*m0/m
  V_Z <-  var(Z1)*m1/m + var(Z0)*m0/m
  Cov_YZ <- cov(Y1,Z1)*m1/m + cov(Y0,Z0)*m0/m
  variance <- m^3/m1/m0 * (V_Y/n_Y^2 + V_Z/n_Z^2 - 2 * Cov_YZ/n_Y/n_Z)
  # variance <- m^3/m1/m0 * (V_Y/n_Y^2 + V_Z/n_Z^2) - 
  #   2 * Cov_YZ * m1 * (m-m1)/m * n_Y * n_Z / sum(Y1) / sum(Y0) / sum(Z1) / sum(Z0)
  CI.lower <- qnorm(0.025, mean = est, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = est, sd = sqrt(variance))
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sqrt(variance), CI.lower, CI.upper, power, coverage))
}

test_pos_frac_est <- function(Y, Z, A, lambda = 1){
  m <- length(Y)
  m1 <- sum(A)
  m0 <- m - m1
  a <- Y/(Y+Z)
  r <- sum(Z)/sum(Y)
  TT <- mean(a[A==1]) - mean(a[A==0])
  U <- Z/Y * (A*lambda + 1-A)
  true_mean <- mean((lambda-1)*U/(lambda+U)/(1+U))
  target_mean <- 2 * r * (lambda^2 - 1)/((2+r)*lambda + r)/(r*lambda + 2 + r)
  variance <- var(a[A==1])/m1 + var(a[A==0])/m0
  f <- function(l) {2 * r * (l^2 - 1) - TT * ((2+r)*l + r) * (r*l + 2 + r)}
  est <- tryCatch(uniroot(f, c(0,2), tol = 1e-6)$root,
                  error = function(e){NA})
  CI.lower <- qnorm(0.025, mean = TT, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = TT, sd = sqrt(variance))
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= target_mean & CI.upper >= target_mean
  return(c(log(est)-log(lambda), target_mean-true_mean, CI.lower, CI.upper, power, coverage))
}


log_contrast_est <- function(Y, Z, A, lambda = 1){
  m <- length(Y)
  m1 <- sum(A)
  m0 <- m - m1
  if(min(Y) == 0) {return(NA)}
  L <- log(Y/Z)
  est <- mean(L[A==1]) - mean(L[A==0])
  variance <- var(L[A==1])/m1 + var(L[A==0])/m0
  CI.lower <- qnorm(0.025, mean = est, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = est, sd = sqrt(variance))
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sqrt(variance), CI.lower, CI.upper, power, coverage))
}

cov_adj_est <- function(Y, Z, A, X = NULL, lambda = 1){
  m <- length(Y)
  m1 <- sum(A)
  m0 <- m - m1
  if(min(Y) == 0) {return(NA)}
  L <- log(Y/Z)
  X <- as.matrix(X)
  V_X <- cov(X)
  beta <- solve(V_X) %*% (cov(X[A==1,], L[A==1])*m1/m + cov(X[A==0,], L[A==0])*m0/m)
  est <-  mean(L[A==1]) - mean(L[A==0]) - t(beta) %*% (colMeans(X[A==1,,drop=F]) - colMeans(X[A==0,,drop=F]))
  variance <- var(L[A==1] - X[A==1,,drop=F] %*% beta)/m1 * (m1-1)/(m1-1-ncol(X)) +
    var(L[A==0] - X[A==0,,drop=F] %*% beta)/m0 * (m0-1)/(m0-1-ncol(X))
  CI.lower <- qnorm(0.025, mean = est, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = est, sd = sqrt(variance))
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sqrt(variance), CI.lower, CI.upper, power, coverage))
}

mixed_est <- function(Y, Z, A, X = NULL, lambda = 1){
  d <- data.frame(Y = Y, Z = Z, A = A, cluster = as.factor(1:length(Y)))
  if(!is.null(X)) {
    d <- cbind(d, X)
    formu <- as.formula(paste("cbind(Y, Z) ~ A + (1|cluster)",paste(colnames(X), collapse = "+"), sep="+"))
  }else{
    formu <- as.formula("cbind(Y, Z) ~ A + (1|cluster)")
  }
  temp_fit <- summary(glmer(formu, data = d, family = binomial))$coefficients
  est <- temp_fit[2,1]
  sd <- temp_fit[2,2]
  CI.lower <- qnorm(0.025, mean = est, sd = sd)
  CI.upper <- qnorm(0.975, mean = est, sd = sd)
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sd, CI.lower, CI.upper, reject_null, coverage))
}
gee_est <- function(Y, Z, A, X = NULL, lambda = 1){
  d <- data.frame(Y = Y, Z = Z, A = A, cluster = as.factor(1:length(Y)))
  if(!is.null(X)) {d <- cbind(d, X)}
  gee_fit <- geeglm(cbind(Y, Z) ~ .-cluster, data = d, id = cluster, family = binomial, corstr = "exchangeable")
  gee_coef <- summary(gee_fit)$coefficients
  est <- gee_coef[2,1]
  sd <- gee_coef[2,2]
  CI.lower <- qnorm(0.025, mean = est, sd = sd)
  CI.upper <- qnorm(0.975, mean = est, sd = sd)
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sd, CI.lower, CI.upper, reject_null, coverage))
}


# data analysis -------------
Y <- d$test_pos
Z <- d$test_neg
A <- d$intervention
odds_ratio_est(Y, Z, A, 1) %>% exp
log_contrast_est(Y, Z, A, 1) %>% exp
test_pos_frac_est(Y, Z, A, 1)
cov_adj_est(Y,Z,A, X = cbind(X1,X2,population, prop_children),1) %>% exp
mixed_est(Y,Z,A, X = scale(cbind(X1,X2,population, prop_children)), 1) %>% exp
gee_est(Y,Z,A, X= cbind(X1,X2,population, prop_children), 1) %>% exp

ll <- seq(exp(-5),exp(1), by = 0.01)
xx <- as.matrix( map_dfc(ll, ~test_pos_frac_est(Y, Z, A, .)))
range(ll[which(xx[6,]==1)])
# dose-response relationship--------------
beta_grid <- seq(-5, 3, by = 0.02)
L <- log(d$test_pos/d$test_neg)
D <- d$wei
A <- d$intervention
m <- length(A)
m1 <- sum(A)
p_value <- map_dbl(beta_grid, function(beta0){
  est <- mean(L[A==1] - beta0 * D[A==1]) - mean(L[A==0] - beta0 * D[A==0])
  var <- var(L[A==1] - beta0 * D[A==1])/m1 + var(L[A==0] - beta0 * D[A==0])/(m - m1)
  2 * pt(-abs(est)/sqrt(var), df = m)
})

range(beta_grid[p_value >= 0.05])
beta_grid[which.max(p_value)]