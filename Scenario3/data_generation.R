expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

backtrack <- function(x, dx, f, Df, al=0.2, be=0.6) {
  t <- 1
  g <- Df(x)
  u <- al*(t(g)%*%dx)
  k <- 1
  repeat {
    if (f(x + t*dx) <= f(x) + t*u) break
    t <- be*t
    # print(paste0("backtrack ",k))
    k <- k + 1
  }
  return(t)
}

update_laplacePoisson <- function(Y_laplacePoisson, X_laplacePoisson,alpha){
  iter = 0
  max.iter = 1000
  eps = 0.0001
  made.changes = TRUE
  len <- length(Y_laplacePoisson)
  # fit <- glm(Y_laplacePoisson ~ -1+X_laplacePoisson,family = poisson("log"))
  mu_laplacePoisson <- rnorm(4, 0, 1)
  #log likelihood
  f = function(b){-sum(Y_laplacePoisson*(X_laplacePoisson%*%b) - exp(X_laplacePoisson%*%b))} #negative log(posterior) ~ log(normal prior*likelihood)
  Df = function(b){-t(X_laplacePoisson)%*%(Y_laplacePoisson - exp(X_laplacePoisson%*%b))} #first derivative
  
  while(made.changes & (iter < max.iter)) {
    iter = iter + 1
    made.changes = FALSE
    H <- (t(X_laplacePoisson)%*%diag(c(exp(X_laplacePoisson%*%mu_laplacePoisson)), nrow = len)%*%X_laplacePoisson)#heissan matrix
    
    d <- -ginv(H)%*%Df(mu_laplacePoisson) #direction
    a = backtrack(mu_laplacePoisson, d, f, Df, al=0.2, be=0.8)
    new.beta_t <- mu_laplacePoisson
    mu_laplacePoisson <- mu_laplacePoisson + a*d
    relative.change = max(abs(new.beta_t - mu_laplacePoisson))
    made.changes = (relative.change > eps)
  }
  
  if (made.changes) {
    warning("Newton's method terminated before convergence")
  }
  H <- (t(X_laplacePoisson)%*%diag(c(exp(X_laplacePoisson%*%mu_laplacePoisson)), nrow = len)%*%X_laplacePoisson)
  variance <- ginv(H)
  # beta_laplacePoisson <- mvrnorm(1, mu_laplacePoisson, variance)
  return(list("mu_laplacePoisson" = mu_laplacePoisson, "variance" = variance))
}

apply_TS <- function(dta, t){
  prob_a <- c()
  for (i in (1:sample_size)) {
    data <- dta[dta$userid == i,]
    Y_use <- data$Y[1:(t-1)]
    A_use <- data$A[1:(t-1)]
    S_use <- data$S[1:(t-1)]
    X_use <- cbind(rep(1,t-1), A_use, S_use, A_use * S_use)
    X_com1 <- matrix(c(1, 1, data$S[t], data$S[t]), nrow = 1, ncol = 4)
    X_com0 <- matrix(c(1, 0, data$S[t], 0), nrow = 1, ncol = 4)
    alpha <- 1
    result <- update_laplacePoisson(Y_use, X_use, alpha)
    beta <- result$mu_laplacePoisson
    variance <- result$variance
    prob_a[i] <- max(0.05, min(0.95, pnorm(0, mean = (X_com0 - X_com1) %*% beta, 
                       sd = (X_com0 - X_com1) %*% variance %*% t(X_com0 - X_com1))))
  }
  return(prob_a)
}

generate_data <- function(sample_size, total_T) {
  
  baseline_Y_S0 <- 2.2
  baseline_Y_S1 <- 2.5
  baseline_Y_S2 <- 2.4
  
  beta_0 <- 0.1
  beta_1 <- 0.3
  
  eta1 <- -0.5
  eta2 <- 0.5
  
  gam1 <- -0.4
  gam2 <- 0.1
  
  
  df_names <- c("userid", "day", "Y", "A", "S", "S2", "prob_zero","prob_Y", "prob_Y_A0", "prob_A")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  
  A <- rep(0, sample_size)
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0)
    if (t <= 20) {
      prob_a <- rep(0.5, sample_size)
    } else {
      prob_a <- apply_TS(dta, t)
    }
    dta$prob_A[row_index] <- prob_a
    A <- rbinom(sample_size, 1, dta$prob_A[row_index])
    dta$A[row_index] <- A
    
    
        ### ZINB ###
        dta$prob_Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                           ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
        dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * (beta_0 + beta_1 * dta$S[row_index]))
        dta$prob_zero[row_index] <- exp(gam1* (dta$S[row_index]+0.1) + gam2 * dta$S[row_index] * dta$A[row_index])
        
        zero <- rbinom(sample_size, 1, dta$prob_zero[row_index])
        dta$Y[row_index] <- zero * rnbinom(sample_size, size = 1, mu = dta$prob_Y[row_index])
    
  }
  
  return(dta)
}

## true beta for (Intercept, S)
beta_true <- c(0.1, 0.4)

## true beta for Intercept
beta_true_marginal <- 0.4599



# compute marginal beta_true
if (0) {

  ### analytically, we have ###
  
  baseline_Y_S0 <- 2.2
  baseline_Y_S1 <- 2.5
  baseline_Y_S2 <- 2.4
  
  beta_0 <- 0.1
  beta_1 <- 0.3
  
  # prob_a <- 0.2
  
  gam1 <- -0.4
  gam2 <- 0.1
  
  numerator <- exp(gam1 * 0.1) * baseline_Y_S0 * exp(beta_0) + exp(gam1 * 1.1 + gam2) * baseline_Y_S1 * exp(beta_0 + beta_1) + exp(gam1 * 2.1 + 2 * gam2) * baseline_Y_S2 * exp(beta_0 + 2 * beta_1)
  denominator <- exp(gam1 * 0.1) * baseline_Y_S0 + exp(gam1 * 1.1) * baseline_Y_S1 + exp(gam1 * 2.1) * baseline_Y_S2
  log(numerator / denominator)

}


