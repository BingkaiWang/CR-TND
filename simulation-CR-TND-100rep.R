rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack)
library(lme4)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

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

cov_adj_est <- function(Y, Z, A, X, lambda = 1){
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
lambda <- 0.6
sim_size <- 10000

s_result <- foreach(i = 1:100, .packages = c("tidyverse", "lme4", "geepack")) %dopar% {
  alpha <- rbeta(n = 24, shape1 = 0.5, shape2 = 0.5)
  simulated_data_CR <- map(1:sim_size, function(j){
    sim_Y0 <- rmultinom(1, size = sum(Y0), prob = Y0/sum(Y0)) * (2 * population)
    sim_Z0 <- rmultinom(1, size = sum(Z0), prob = Z0/sum(Z0)) / (2 * population)
    sim_Z1 <- alpha * sim_Z0
    sim_Y1 <- lambda * alpha * sim_Y0
    A <- sample(rep(0:1, each = 12), size = 24, replace = F)
    data.frame(A = A,
               Y = A * sim_Y1 + (1-A) * sim_Y0,
               Z = A * sim_Z1 + (1-A) * sim_Z0)
  })
  
    x1 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~odds_ratio_est(.$Y, .$Z, .$A, lambda = lambda))))
    x2 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~log_contrast_est(.$Y, .$Z, .$A, lambda = lambda))))
    x3 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~cov_adj_est(.$Y, .$Z, .$A, cbind(population), lambda = lambda))))
    x4 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~test_pos_frac_est(.$Y, .$Z, .$A, lambda = lambda))))
    x5 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~mixed_est(.$Y, .$Z, .$A, cbind(population), lambda = lambda))))
    x6 <- suppressMessages(as.matrix(map_dfc(simulated_data_CR, ~gee_est(.$Y, .$Z, .$A,cbind(population), lambda = lambda))))

  summary <- data.frame(bias = map_dbl(list(x1,x2,x3,x4,x5,x6), ~mean(.[1,], na.rm=T)),
                        se = map_dbl(list(x1,x2,x3,x4,x5,x6), ~sd(.[1,], na.rm=T)),
                        ese = map_dbl(list(x1,x2,x3,x4,x5,x6), ~mean(.[2,], na.rm=T)),
                        power = map_dbl(list(x1,x2,x3,x4,x5,x6), ~mean(.[5,], na.rm=T)),
                        CP = map_dbl(list(x1,x2,x3,x4,x5,x6), ~mean(.[6,], na.rm=T)))
  rownames(summary) <- c("Odds_ratio", "log_contrast", "cov_adj", "test_pos_frac", "mixed_model", "GEE")
  summary
}

stopCluster(cl)

# visualization ----------------
library(ggpattern)
bias_summary <- t(as.matrix(map_dfc(s_result, ~.[,1])))
colnames(bias_summary) <- c("Odds\ ratio", "Log-contrast", "Covariate-adjusted", "Test-positive\ fraction", "GLMM", "GEE")
bias_summary <- abs(bias_summary) %>% data.frame(check.names = FALSE) %>%
  pivot_longer(`Odds ratio`:GEE, names_to = "Estimator", values_to = "Absolute_bias") %>%
  mutate(bin = cut_width(Absolute_bias, width = 0.04, center = 0.02)) %>%
  mutate(Estimator = factor(Estimator, levels = c("Odds ratio", "Test-positive fraction", "GLMM", "GEE", "Log-contrast", "Covariate-adjusted")))
levels(bias_summary$bin) <- c(levels(bias_summary$bin)[1:5], rep("> 0.2",6))
p_bias <- ggplot(bias_summary) + 
  geom_bar_pattern(aes(bin, pattern = Estimator, pattern_angle = Estimator, pattern_fill = Estimator,  pattern_spacing = Estimator),
                   colour = 'black', fill = 'white', width= 0.6 ,stat="count",position = "dodge") + 
  scale_pattern_fill_manual(values = c("black","white","black","white", "black","black")) +
  scale_pattern_spacing_manual(values = c(0.02,0.05,0.02,0.05,0.01, 0.02)) +
  theme_bw() +
  theme(text = element_text(size =30, family = "Times"), legend.position = "bottom", 
        legend.background = element_rect(color = "black"), legend.text = element_text(size = 30), legend.key.size = unit(1, "cm")) +
  labs(x = "Absolute bias", y = "Percentage (%)")

cp_summary <- t(as.matrix(map_dfc(s_result, ~.[,5])))
colnames(cp_summary) <- c("Odds\ ratio", "Log-contrast", "Covariate-adjusted", "Test-positive\ fraction", "GLMM", "GEE")
cp_summary <- abs(cp_summary) %>% data.frame(check.names = FALSE) %>%
  pivot_longer(`Odds ratio`:GEE, names_to = "Estimator", values_to = "CP") %>%
  mutate(bin = cut_width(CP, width = 0.04, center = 0.95)) %>%
  mutate(Estimator = factor(Estimator, levels = c("Odds ratio", "Test-positive fraction", "GLMM", "GEE", "Log-contrast", "Covariate-adjusted")))
levels(cp_summary$bin) <- c(levels(cp_summary$bin)[1:5], levels(cp_summary$bin)[5])
p_cp <- ggplot(cp_summary) + 
  geom_bar_pattern(aes(bin, pattern = Estimator, pattern_angle = Estimator, pattern_fill = Estimator,  pattern_spacing = Estimator), 
                   colour = 'black', fill = 'white', width= 0.6 ,stat="count",position = "dodge") + 
  scale_pattern_fill_manual(values = c("black","white","black","white", "black","black")) +
  scale_pattern_spacing_manual(values = c(0.02,0.05,0.02,0.05,0.01, 0.02)) +
  theme_bw() +
  theme(text = element_text(size =30, family = "Times"), legend.position = "bottom", legend.background = element_rect(color = "black")) +
  labs(x = "Coverage Probability", y = "Percentage (%)")
# p_combined <- cowplot::plot_grid(p_bias, p_cp, ncol = 2)
# cowplot::save_plot(plot = p_combined, filename = "figure1-blackwhite.png", base_height = 10, base_asp = 2)
leg <- cowplot::get_legend(p_bias)
p_combined <- cowplot::plot_grid(p_bias+theme(legend.position='none'), p_cp+theme(legend.position='none'), ncol = 2)
pp <- cowplot::plot_grid(p_combined, leg, ncol = 1, rel_heights = c(1, 0.12))
cowplot::save_plot(plot = pp, filename = "figure1-blackwhite.png", base_height = 10, base_asp = 2.2)

