rm(list = ls())
set.seed(1234)
library(tidyverse)
library(staggered)
library(geepack)
library(lme4)
library(tictoc)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

# estimators ------
oddsratio_SWD <- function(d, lambda = 1, n_fisher_permutation = 500){
  A <- d[,1]
  dengue_data1 <- dengue_data0 <- d[ ,2:10]
  OFI_data1 <- OFI_data0 <- d[ ,11:19]
  if(min(dengue_data1) == 0 | min(OFI_data1) == 0){
    return(rep(NA, 6))
  }
  for(i in 1:length(A)){
    dengue_data1[i, 1:(A[i]-1)] <- NA
    dengue_data0[i, A[i]:9] <- NA
    OFI_data1[i, 1:(A[i]-1)] <- NA
    OFI_data0[i, A[i]:9] <- NA
  }
  dengue_data1 <- as.vector(dengue_data1[,2:8]) %>% .[!is.na(.)]
  dengue_data0 <- as.vector(dengue_data0[,2:8]) %>% .[!is.na(.)]
  OFI_data1 <- as.vector(OFI_data1[,2:8]) %>% .[!is.na(.)]
  OFI_data0 <- as.vector(OFI_data0[,2:8]) %>% .[!is.na(.)]
  est <- log(sum(dengue_data1)/sum(dengue_data0)) - log(sum(OFI_data1)/sum(OFI_data0))
  est_perm <- map_dbl(1:n_fisher_permutation, function(j){
    AA <- sample(A, size = length(A))
    dengue_data1 <- dengue_data0 <- d[ ,2:10]
    OFI_data1 <- OFI_data0 <- d[ ,11:19]
    for(i in 1:length(A)){
      dengue_data1[i, 1:(AA[i]-1)] <- NA
      dengue_data0[i, AA[i]:9] <- NA
      OFI_data1[i, 1:(AA[i]-1)] <- NA
      OFI_data0[i, AA[i]:9] <- NA
    }
    dengue_data1 <- as.vector(dengue_data1[,2:8]) %>% .[!is.na(.)]
    dengue_data0 <- as.vector(dengue_data0[,2:8]) %>% .[!is.na(.)]
    OFI_data1 <- as.vector(OFI_data1[,2:8]) %>% .[!is.na(.)]
    OFI_data0 <- as.vector(OFI_data0[,2:8]) %>% .[!is.na(.)]
    log(sum(dengue_data1)/sum(dengue_data0)) - log(sum(OFI_data1)/sum(OFI_data0))
  })
  CI.lower <- qnorm(0.025, mean = est, sd = sd(est_perm))
  CI.upper <- qnorm(0.975, mean = est, sd = sd(est_perm))
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sd(est_perm), CI.lower, CI.upper, reject_null, coverage))
}
est2_SWD <- function(d, lambda = 1){
  A <- d[,1]
  m <- map_dbl(1:7, ~sum(A<=.+1))
  dengue_data <- dengue_data0 <- d[ ,3:9]
  OFI_data <- OFI_data0 <- d[ ,12:18]
  if(min(dengue_data) == 0 | min(OFI_data) == 0){
    return(rep(NA, 6))
  }
  Y <- log(dengue_data/OFI_data)
  S_mat <- cov(Y, use = "all.obs")
  for(i in 1:ncol(S_mat)){
    for(j in 1:ncol(S_mat)){
      S_mat[i,j] <- nrow(Y)/m[max(i,j)]/(nrow(Y)-m[min(i,j)]) * S_mat[i,j]
    }
  }
  if(min(eigen(S_mat)$values) == 0){
    weights <- rep(1/7, 7)
    variance <- sum(S_mat)/49
  }else{
    S_inv <- solve(S_mat)
    weights <- rowSums(S_inv)/sum(S_inv)
    variance <- 1/sum(S_inv)
  }
  est <- map_dbl(1:7, function(j){
    sum(Y[A <= j+1, j])/m[j] - sum(Y[A > j+1, j])/(length(A)-m[j])
  }) %*% weights
  CI.lower <- qnorm(0.025, mean = est, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = est, sd = sqrt(variance))
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sqrt(variance), CI.lower, CI.upper, reject_null, coverage))
}
swd_est <- function(Y, Z, A, lambda = 1, optim = T){
  TT <- ncol(Y)
  m <- nrow(Y)
  A_mat <- matrix(0, nrow = m, ncol = TT)
  for(i in 1:m){
    A_mat[i, A[i]: TT] <- 1
  }
  m_t <- colSums(A_mat == 1)
  if(min(Y) == 0) {return(NA)}
  L <- log(Y/Z)
  if(optim == T){
    L0 <- L - log(lambda) * A_mat
    S <- cov(L0, use = "all.obs")[1:(TT-1), 1:(TT-1)]
    for(i in 1:ncol(S)){
      for(j in 1:ncol(S)){
        S[i,j] <- m/m_t[max(i,j)]/(m-m_t[min(i,j)]) * S[i,j]
      }
    }
    if(min(eigen(S)$values) <= 1e-5){
      weights <- rep(1/(TT-1), TT-1)
      variance <- sum(S)/(TT-1)^2
    }else{
      S_inv <- solve(S)
      weights <- rowSums(S_inv)/sum(S_inv)
      variance <- 1/sum(S_inv)
    }
  } else {
    S <- matrix(NA, nrow = TT-1, ncol = TT-1)
    for(t1 in 1:(TT-1)){
      for(t2 in t1:(TT-1)){
        index11 <- which(A_mat[,t1] == 1)
        index01 <- which(A_mat[,t2] == 1 & A_mat[,t1] == 0)
        index00 <- which(A_mat[,t2] == 0)
        if(0){
          S[t1, t2] <- cov(L[index11, t1], L[index11, t2])/m_t[t2] +
            cov(L[index00, t1], L[index00, t2])/(m-m_t[t1])
        } else{
          if(length(index11) >= length(index01) & length(index11) >= length(index00)){
            S[t1, t2] <- S[t2, t1] <- m/m_t[t2]/(m-m_t[t1]) * cov(L[index11, t1], L[index11, t2])
          } else if(length(index01) >= length(index11) & length(index01) >= length(index00)) {
            S[t1, t2] <- S[t2, t1] <- m/m_t[t2]/(m-m_t[t1]) * cov(L[index01, t1], L[index01, t2])
          } else {
            S[t1, t2] <- S[t2, t1] <- m/m_t[t2]/(m-m_t[t1]) * cov(L[index00, t1], L[index00, t2])
          }
        }
        
      }
    }
    weights <- rep(1/(TT-1), TT-1)
    variance <- sum(S)/(TT-1)^2
  }
  est <- map_dbl(1:(TT-1), function(j){
    mean(L[A <= j, j]) - mean(L[A > j, j])
  }) %*% weights
  CI.lower <- qnorm(0.025, mean = est, sd = sqrt(variance))
  CI.upper <- qnorm(0.975, mean = est, sd = sqrt(variance))
  power <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sqrt(variance), CI.lower, CI.upper, power, coverage))
}
mixed_est_SWD <- function(d, lambda = 1){
  transformed_data <- pivot_longer(mutate(data.frame(d), cluster = 1:24), Y03dengue:Y13OFI, names_to = "time", values_to = "Outcomes")
  time_index <- c(3, 5, 6, 7, 8, 10, 11, 12, 13)
  transformed_data <- transformed_data %>%
    mutate(disease = ifelse(str_detect(transformed_data$time, "dengue"), "dengue", "OFI")) %>%
    mutate(time = as.numeric(str_extract(transformed_data$time, "[0-9]+"))) %>%
    mutate(treatment = time_index[transformed_data$A] <= time) %>%
    filter(!(time %in% c(3, 13))) %>%
    mutate(time = as.factor(time)) %>%
    mutate(cluster = as.factor(cluster)) %>%
    pivot_wider(names_from = disease, values_from = Outcomes)
  temp_fit <- glmer(cbind(dengue, OFI) ~ treatment + time + (1|cluster) + (1|cluster:time), data = transformed_data, family = binomial)
  temp_coef <- summary(temp_fit)$coef
  est <- temp_coef[2,1]
  sd <- temp_coef[2,2]
  CI.lower <- qnorm(0.025, mean = est, sd = sd)
  CI.upper <- qnorm(0.975, mean = est, sd = sd)
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est - log(lambda), sd, CI.lower, CI.upper, reject_null, coverage))
}
gee_est_SWD <- function(d, lambda = 1){
  transformed_data <- pivot_longer(mutate(data.frame(d), cluster = 1:24), Y03dengue:Y13OFI, names_to = "time", values_to = "Outcomes")
  time_index <- c(3, 5, 6, 7, 8, 10, 11, 12, 13)
  transformed_data <- transformed_data %>%
    mutate(disease = ifelse(str_detect(transformed_data$time, "dengue"), "dengue", "OFI")) %>%
    mutate(time = as.numeric(str_extract(transformed_data$time, "[0-9]+"))) %>%
    mutate(treatment = time_index[transformed_data$A] <= time) %>%
    mutate(time = as.factor(time)) %>%
    mutate(cluster = as.factor(cluster)) %>%
    pivot_wider(names_from = disease, values_from = Outcomes)  
  temp_fit <- geeglm(cbind(dengue, OFI) ~ treatment, data = transformed_data, id = cluster, family = binomial, corstr = "exchangeable")
  temp_coef <- summary(temp_fit)$coef
  est <- temp_coef[2,1]
  sd <- temp_coef[2,2]
  CI.lower <- qnorm(0.025, mean = est, sd = sd)
  CI.upper <- qnorm(0.975, mean = est, sd = sd)
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sd, CI.lower, CI.upper, reject_null, coverage))
}
staggered_est <- function(d){
  transformed_data <- pivot_longer(mutate(data.frame(d), cluster = 1:24), Y03dengue:Y13OFI, names_to = "time", values_to = "Outcomes")
  time_index <- c(3, 5, 6, 7, 8, 10, 11, 12, 13)
  if(min(transformed_data$Outcomes) == 0){return(NA)}
  transformed_data <- transformed_data %>%
    mutate(disease = ifelse(str_detect(transformed_data$time, "dengue"), "dengue", "OFI")) %>%
    mutate(time = as.numeric(str_extract(transformed_data$time, "[0-9]+"))) %>%
    mutate(treatment = time_index[transformed_data$A] <= time) %>%
    mutate(time = as.numeric(as.factor(time))) %>%
    mutate(cluster = as.factor(cluster)) %>%
    pivot_wider(names_from = disease, values_from = Outcomes) %>%
    mutate(L = log(dengue/OFI))
  result <- staggered(df = transformed_data, 
                      i = "cluster",
                      t = "time",
                      g = "A",
                      y = "L", beta = 0,
                      estimand = "simple") %>% t
  est <- result[1]
  sd <- result[2]
  CI.lower <- qnorm(0.025, mean = est, sd = sd)
  CI.upper <- qnorm(0.975, mean = est, sd = sd)
  reject_null <- CI.lower > 0 | CI.upper < 0
  coverage <- CI.lower <= log(lambda) & CI.upper >= log(lambda)
  return(c(est-log(lambda), sd, CI.lower, CI.upper, reject_null, coverage))
}


# simulation -------
dengue_cases <- c(
  1, 13, 19, 37, 29, 42, 48, 18, 26, 34,
  2, 14, 14, 30, 27, 34, 37, 15, 25, 34,
  3, 35, 32, 39, 43, 62, 52, 25, 40, 38,
  4, 9, 13, 13, 8, 18, 18, 6, 7, 9,
  5, 17, 25, 69, 60, 36, 53, 34, 47, 71,
  6, 37, 38, 77, 72, 75, 89, 84, 120, 104,
  7, 23, 28, 48, 51, 85, 76, 28, 40, 36,
  8, 20, 32, 51, 57, 66, 41, 13, 36, 37,
  9, 25, 29, 46, 41, 57, 48, 15, 27, 25,
  10, 14, 25, 53, 49, 41, 31, 9, 35, 42,
  11, 40, 61, 78, 64, 84, 98, 57, 62, 71,
  12, 33, 54, 74, 59, 80, 80, 44, 63, 69,
  13, 35, 52, 79, 86, 119, 112, 49, 56, 76,
  14, 28, 39, 57, 48, 59, 56, 29, 49, 62,
  15, 30, 39, 56, 46, 52, 40, 20, 25, 27,
  16, 22, 51, 68, 47, 56, 43, 19, 36, 38,
  17, 12, 18, 25, 22, 20, 14, 8, 17, 16,
  18, 41, 55, 112, 93, 130, 151, 81, 139, 128,
  19, 16, 27, 69, 71, 53, 44, 24, 47, 69,
  20, 19, 37, 43, 28, 45, 41, 30, 79, 77,
  21, 24, 45, 63, 49, 59, 62, 42, 73, 68,
  22, 33, 57, 72, 59, 84, 73, 35, 66, 62,
  23, 12, 19, 29, 29, 36, 29, 14, 34, 32,
  24, 21, 40, 67, 90, 151, 106, 27, 72, 76)
OFI_cases14 <- c(
  1, 486,
  2, 155,
  3, 1197,
  4, 255,
  5, 249,
  6, 710,
  7, 658,
  8, 714,
  9, 478,
  10, 376,
  11, 388,
  12, 426,
  13, 842,
  14, 547,
  15, 285,
  16, 586,
  17, 344,
  18, 484,
  19, 151,
  20, 223,
  21, 522,
  22, 804,
  23, 286,
  24, 792
)
dengue_cases <- t(matrix(dengue_cases, ncol = 24)) %>% data.frame
colnames(dengue_cases) <- c("cluster_id", "Y03", "Y05", "Y06", "Y07", "Y08", "Y10", "Y11", "Y12", "Y13")
OFI_cases14 <-  t(matrix(OFI_cases14, ncol = 24)) %>% data.frame
colnames(OFI_cases14) <- c("cluster_id", "OFI14")
sim_size <- 10000
alpha <- rbeta(n = 9*24, shape1 = 0.5, shape2 = 0.5) %>% matrix(nrow = 24, ncol = 9)
lambda <- 0.2
ksi <- eta <- 0.1 # rates of imprecise tests
n_dengue <- colSums(dengue_cases[,-1])
dengue_p <- apply(dengue_cases[,-1], 2, function(x){x/sum(x)})
n_OFI <- round(rnorm(sum(OFI_cases14[,-1])/n_dengue[9],mean=10,sd = 1) * n_dengue)
OFI_p <- OFI_cases14[,-1]/sum(OFI_cases14[,-1])
simulated_data_SW <- map(1:sim_size, function(j){
  sim_OFI_cases0 <- as.matrix(map_dfc(n_OFI, ~rmultinom(1, size = ., prob = OFI_cases14[,-1]/sum(OFI_cases14[,-1]))))
  sim_dengue_cases0 <- apply(dengue_cases[,-1], 2, function(x){rmultinom(1,size=sum(x), prob = x/sum(x))})
  sim_OFI_cases1 <- alpha * sim_OFI_cases0 
  sim_dengue_cases1 <- lambda * alpha * sim_dengue_cases0
  obs_dengue_cases0 <- (1-ksi) * sim_dengue_cases0 + eta * sim_OFI_cases0
  obs_dengue_cases1 <- (1-ksi) * sim_dengue_cases1 + eta * sim_OFI_cases1
  obs_OFI_cases0 <- ksi * sim_dengue_cases0 + (1-eta) * sim_OFI_cases0
  obs_OFI_cases1 <- ksi * sim_dengue_cases1 + (1-eta) * sim_OFI_cases1
  A <- sample(rep(2:9, each = 3), size = 24, replace = F)
  sim_OFI <- as.matrix(map_dfr(1:24, function(i){
    temp <- obs_OFI_cases0[i,]
    temp[A[i]:9] <- obs_OFI_cases1[i,A[i]:9]
    return(temp)
  }))
  colnames(sim_OFI) <- paste0(colnames(sim_OFI), "OFI")
  sim_dengue <- as.matrix(map_dfr(1:24, function(i){
    temp <- obs_dengue_cases0[i,]
    temp[A[i]:9] <- obs_dengue_cases1[i,A[i]:9]
    return(temp)
  }))
  colnames(sim_dengue) <- paste0(colnames(sim_dengue), "dengue")
  cbind(A, sim_dengue, sim_OFI)
})


log_contrast_result_optim <- foreach(iter = 1:sim_size, .combine = cbind, .packages = c("tidyverse")) %dopar% {
  A <- simulated_data_SW[[iter]][,1] - 1
  Y <- simulated_data_SW[[iter]][,3:10]
  Z <- simulated_data_SW[[iter]][,12:19]
  swd_est(Y,Z,A, lambda = lambda, optim = T)
}
log_contrast_result <- foreach(iter = 1:sim_size, .combine = cbind, .packages = c("tidyverse")) %dopar% {
  A <- simulated_data_SW[[iter]][,1] - 1
  Y <- simulated_data_SW[[iter]][,3:10]
  Z <- simulated_data_SW[[iter]][,12:19]
  swd_est(Y,Z,A, lambda = lambda, optim = F)
}
mm_result <- foreach(i = 1:sim_size, .combine = cbind, .packages = c("tidyverse", "lme4")) %dopar% {
  mixed_est_SWD(simulated_data_SW[[i]], lambda = lambda)
}
gee_result <- foreach(i = 1:sim_size, .combine = cbind, .packages = c("tidyverse", "geepack")) %dopar% {
  gee_est_SWD(simulated_data_SW[[i]], lambda = lambda)
}
stopCluster(cl)


result <- list(log_contrast_result_optim, log_contrast_result, mm_result, gee_result)
summary <- data.frame(bias = map_dbl(result, ~mean(.[1,], na.rm=T)),
                      se = map_dbl(result, ~sd(.[1,], na.rm=T)),
                      ese = map_dbl(result, ~mean(.[2,], na.rm=T)),
                      power = map_dbl(result, ~mean(.[5,], na.rm=T)),
                      CP = map_dbl(result, ~mean(.[6,], na.rm=T)))
rownames(summary) <- c("log_contrast_optim", "log_contrast", "mixed", "GEE")
xtable::xtable(summary)