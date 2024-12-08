library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package

get_alpha_beta_from_multiroot_result <- function(root, p, q, num_trt)
{
  if (p == 1) {
    beta_root <- root$root[(q+1):(q+num_trt)]
    beta_root <- t(as.matrix(beta_root))
  } else {
    beta_root <- matrix(NA, nrow = p, ncol = num_trt)
    for (k in 1:num_trt) {
      beta_root[,k] <- root$root[(q + 1 + (k-1)*p) : (q + k*p)]
    }
  }
  if (q == 1) {
    alpha_root <- root$root[1]
  } else {
    alpha_root <- as.matrix(root$root[1:q])
    #alpha_root <- as.vector(root$root[1:q])
  }
  return(list(alpha = alpha_root, beta = beta_root))
}

find_change_location <- function(v){
  n <- length(v)
  if (n <= 1) {
    stop("The vector need to have length > 1.")
  }
  return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}



fit_EMEE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
    estimator_initial_value = NULL
)
{


  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  A <- dta[, treatment_varname]
  # checking for NA in treatment indicator    
  if (any(is.na(A[avail == 1,]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  A1 <- A[,1]
  A2 <- A[,2]
  
  p_t <- dta[, rand_prob_varname]
  p_t1 <- p_t[,1]
  p_t2 <- p_t[,2]
  cA1 <- A1 - p_t1 # centered A1
  cA2 <- A2 - p_t2 # centered A2

  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde1 <- p_t_tilde2 <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 2) {
      p_t_tilde1 <- rep(rand_prob_tilde[1], total_person_decisionpoint)
      p_t_tilde2 <- rep(rand_prob_tilde[2], total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == 2*total_person_decisionpoint) {
      p_t_tilde1 <- rand_prob_tilde[1:total_person_decisionpoint]
      p_t_tilde2 <- rand_prob_tilde[(total_person_decisionpoint+1):2*total_person_decisionpoint]
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde1 <- A1 - p_t_tilde1
  cA_tilde2 <- A2 - p_t_tilde2
  
  WCLS_weight <- ifelse(A1, p_t_tilde1 / p_t1, ifelse(A2, p_t_tilde2 / p_t2, (1 - p_t_tilde1 - p_t_tilde2) / (1 - p_t1 - p_t2)))
  
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Intercept", control_varname)
  
  ### estimation ###
  
  estimating_equation <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta1 <- as.matrix(theta[(q+1):(q+p)])
    beta2 <- as.matrix(theta[(q + p + 1):(q + 2*p)])
    
    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_AXdm_beta <- exp(A1 * (Xdm %*% beta1) + A2 * (Xdm %*% beta2))
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1)
    
    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * WCLS_weight * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde1 * Xdm[, i])
      ef[q + p + i] <- sum(weight * residual * avail * WCLS_weight * cA_tilde2 * Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = 2 * p + q)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, 2*p + q), msg = cond,
                  f.root = rep(NaN, 2*p + q)))
    })
  
  estimator <- get_alpha_beta_from_multiroot_result(solution, p, q, 2)
  alpha_hat <- as.vector(estimator$alpha)
  beta_hat <- estimator$beta
  
  ### asymptotic variance ###
  
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, 2 * p + q, 2 * p + q))
 
  r_term_collected <- rep(NA, total_person_decisionpoint)
  D_term_collected <- matrix(NA, nrow = 2 * p + q, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = 2 * p + q)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(t(Xdm[it, ]) %*% beta_hat)
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_hat
    } else {
      Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
    }
    
    pre_multiplier <- exp(- A1[it] * Xbeta[1] - A2[it] * Xbeta[2]) * WCLS_weight[it]
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow =  2 * p + q, ncol = 2 * p + q)
    partialD_partialtheta[1:q, 1:q] <- 0
    partialD_partialtheta[1:q, (q + 1):(q + p)] <- -pre_multiplier * A1[it] * (Zdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[1:q, (q + p + 1):(q + 2*p)] <- -pre_multiplier * A2[it] * (Zdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q + 1):(q + 2*p), 1:q] <- 0
    partialD_partialtheta[(q + 1):(q + p), (q + 1):(q + p)] <- -pre_multiplier * A1[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q + 1):(q + p), (q + p + 1):(q + 2*p)] <- -pre_multiplier * A2[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q + p + 1):(q + 2*p), (q + 1):(q + p)] <- -pre_multiplier * A1[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q + p + 1):(q + 2*p), (q + p + 1):(q + 2*p)] <- -pre_multiplier * A2[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    r_term <- (Y[it] - exp(Zalpha + A1[it] * Xbeta[1]+ A2[it] * Xbeta[2])) * avail[it]
    r_term_collected[it] <- r_term
    
    # D_term = D^{(t),T} (dim = (p+q) * 1)
    D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde1[it] * Xdm[it, ], cA_tilde2[it] * Xdm[it, ])
    D_term_collected[, it] <- D_term
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- -exp(Zalpha + A1[it] * Xbeta[1]+ A2[it] * Xbeta[2]) * c(Zdm[it, ], A1[it] * Xdm[it, ], A2[it] * Xdm[it, ]) * avail[it]
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  

  Sigman_summand <- array(NA, dim = c(sample_size, 2 * p + q, 2 * p + q))
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    
    Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  alpha_se <- sqrt(diag(varcov)[1:q])
  beta_se <- sqrt(diag(varcov)[(q + 1):(q + 2 * p)])
  

  
  names(alpha_hat) <- names(alpha_se)  <- Znames
  names(beta_hat) <- names(beta_se) <- rep(Xnames,2)
  
  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}



fit_EMEE_NonP <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
    estimator_initial_value = NULL
)
{

  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  time_length <- length(unique(dta[, decision_time_varname]))
  
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  A <- dta[, treatment_varname]
  # checking for NA in treatment indicator    
  if (any(is.na(A[avail == 1,]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  A1 <- A[,1]
  A2 <- A[,2]
  
  p_t <- dta[, rand_prob_varname]
  p_t1 <- p_t[,1]
  p_t2 <- p_t[,2]
  cA1 <- A1 - p_t1 # centered A1
  cA2 <- A2 - p_t2 # centered A2
  
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde1 <- p_t_tilde2 <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 2) {
      p_t_tilde1 <- rep(rand_prob_tilde[1], total_person_decisionpoint)
      p_t_tilde2 <- rep(rand_prob_tilde[2], total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == 2*total_person_decisionpoint) {
      p_t_tilde1 <- rand_prob_tilde[1:total_person_decisionpoint]
      p_t_tilde2 <- rand_prob_tilde[(total_person_decisionpoint+1):2*total_person_decisionpoint]
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde1 <- A1 - p_t_tilde1
  cA_tilde2 <- A2 - p_t_tilde2
  
  WCLS_weight <- ifelse(A1, p_t_tilde1 / p_t1, ifelse(A2, p_t_tilde2 / p_t2, (1 - p_t_tilde1 - p_t_tilde2) / (1 - p_t1 - p_t2)))
  
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  ### estimation: nuisance parameter [] ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A1==0 & A2 ==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A1==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  df42 = df4[A2==1,]
  df42_prob = df42
  df42_prob$ind = df42$Y > 0
  df42_mean = df42[which(df42_prob$Y > 0),]
  

  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  fit42_prob = gam(formula1, data=df42_prob, family = poisson(link='log'))
  fit42_mean = gam(formula2, data=df42_mean, family = gaussian(link="log"))
  
  
  EY2 = predict(fit42_prob, newdata=df4, type = "response") * predict(fit42_mean, newdata=df4, type = "response")
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  
  ### estimation: parameter of interest [beta] ###
  
  estimating_equation <- function(theta) {
    beta1 <- as.matrix(theta[1:p])
    beta2 <- as.matrix(theta[(p + 1):(2*p)])
    
    exp_AXdm_beta <- exp(A1 * (Xdm %*% beta1) + A2 * (Xdm %*% beta2))
    # exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - p_t_tilde1 - p_t_tilde2) + as.vector(EY1) * p_t_tilde1 * (exp(Xdm %*% beta1))^(-1) + as.vector(EY2) * p_t_tilde2 * (exp(Xdm %*% beta2))^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1)
    
    ef <- rep(NA, 2*p) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde1 * Xdm[, i] )
      ef[p + i] <- sum(weight * residual * avail * WCLS_weight * cA_tilde2 * Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = 2 * p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside estimator 2():")
      message(cond)
      return(list(root = rep(NaN, 2 * p), msg = cond,
                  f.root = rep(NaN, 2 * p)))
    })
  
  beta_hat <- matrix(solution$root, nrow = p)
  
  ### 3. asymptotic variance ###

  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, 2 * p, 2 * p))
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, 2 * p, 2 * p))

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(t(Xdm[it, ]) %*% beta_hat)
    }

    pre_multiplier <- exp(- A1[it] * Xbeta[1] - A2[it] * Xbeta[2]) * WCLS_weight[it]
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow =  2 * p, ncol = 2 * p)
    partialD_partialtheta[1:p, 1:p] <- -pre_multiplier * A1[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ]
    partialD_partialtheta[1:p, (p + 1):(2*p)] <- -pre_multiplier * A2[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(p + 1):(2*p), 1:p] <- -pre_multiplier * A1[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(p + 1):(2*p), (p + 1):(2*p)] <- -pre_multiplier * A2[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - p_t_tilde1[it] - p_t_tilde2[it]) + EY1[it] * p_t_tilde1[it] * exp(-Xbeta[1]) + EY2[it] * p_t_tilde2[it] * exp(-Xbeta[2])
    r_term <- (Y[it] - EH * exp(A1[it] * Xbeta[1]+ A2[it] * Xbeta[2])) * avail[it]
    
    # D_term = D^{(t),T} (dim = (p) * 1)
    D_term <- pre_multiplier * c(cA_tilde1[it] * Xdm[it, ], cA_tilde2[it] * Xdm[it, ])
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- EY1[it] * p_t_tilde1[it] * exp(-Xbeta[1]) * exp(A1[it] * Xbeta[1]+ A2[it] * Xbeta[2]) * avail[it] * Xdm[it, ] +  EY2[it] * p_t_tilde2[it] * exp(-Xbeta[2]) * exp(A1[it] * Xbeta[1]+ A2[it] * Xbeta[2]) * avail[it] * Xdm[it, ]  - EH * exp(A1[it] * Xbeta[1]+ A2[it] * Xbeta[2]) * c(A1[it] * Xdm[it, ], A2[it] * Xdm[it, ]) * avail[it]
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    
    Sigman_summand[it, , ] <- (D_term * r_term) %*% (r_term * t(D_term))
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size / time_length
  Mn_inv <- solve(Mn)
  
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size / time_length
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size / time_length
  beta_se <- sqrt(diag(varcov)[1:(2 * p)])

  
  beta_hat = as.vector(beta_hat)
  names(beta_hat) <- names(beta_se)  <- rep(Xnames, 2)
  
  return(list(beta_hat = beta_hat, beta_se = beta_se,
              varcov = varcov,
              dims = list(p = p),
              f.root = solution$f.root))
}


fit_DR_EMEE_NonP <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
    estimator_initial_value = NULL
)
{

  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  time_length <- length(unique(dta[, decision_time_varname]))
  
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  A <- dta[, treatment_varname]
  # checking for NA in treatment indicator    
  if (any(is.na(A[avail == 1,]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  A1 <- A[,1]
  A2 <- A[,2]
  
  p_t <- dta[, rand_prob_varname]
  p_t1 <- p_t[,1]
  p_t2 <- p_t[,2]
  cA1 <- A1 - p_t1 # centered A1
  cA2 <- A2 - p_t2 # centered A2
  
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde1 <- p_t_tilde2 <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 2) {
      p_t_tilde1 <- rep(rand_prob_tilde[1], total_person_decisionpoint)
      p_t_tilde2 <- rep(rand_prob_tilde[2], total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == 2*total_person_decisionpoint) {
      p_t_tilde1 <- rand_prob_tilde[1:total_person_decisionpoint]
      p_t_tilde2 <- rand_prob_tilde[(total_person_decisionpoint+1):2*total_person_decisionpoint]
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde1 <- A1 - p_t_tilde1
  cA_tilde2 <- A2 - p_t_tilde2
  
  WCLS_weight <- ifelse(A1, p_t_tilde1 / p_t1, ifelse(A2, p_t_tilde2 / p_t2, (1 - p_t_tilde1 - p_t_tilde2) / (1 - p_t1 - p_t2)))
  
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  ### estimation: nuisance parameter [] ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A1==0 & A2 ==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A1==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  df42 = df4[A2==1,]
  df42_prob = df42
  df42_prob$ind = df42$Y > 0
  df42_mean = df42[which(df42_prob$Y > 0),]
  
  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  fit42_prob = gam(formula1, data=df42_prob, family = poisson(link='log'))
  fit42_mean = gam(formula2, data=df42_mean, family = gaussian(link="log"))
  
  
  EY2 = predict(fit42_prob, newdata=df4, type = "response") * predict(fit42_mean, newdata=df4, type = "response")
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  
  ### estimation: parameter of interest [beta] ###
  
  estimating_equation <- function(theta) {
    beta1 <- as.matrix(theta[1:p])
    beta2 <- as.matrix(theta[(p + 1):(2*p)])
    
    exp_AXdm_beta <- exp(A1 * (Xdm %*% beta1) + A2 * (Xdm %*% beta2))
    exp_Xdm_beta1 <- exp(Xdm %*% beta1)
    exp_Xdm_beta2 <- exp(Xdm %*% beta2)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - A1 - A2) + as.vector(EY1) * A1 * (exp(Xdm %*% beta1))^(-1) + as.vector(EY2) * A2 * (exp(Xdm %*% beta2))^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1)
    
    ef <- rep(NA, 2*p) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( (weight * residual * avail * WCLS_weight * cA_tilde1 + 
                       avail * p_t_tilde1 * (1-p_t_tilde1) * (exp_Xdm_beta1^(-1) * as.vector(EY1) - as.vector(EY0)))* Xdm[, i] )
      ef[p + i] <- sum( (weight * residual * avail * WCLS_weight * cA_tilde2 + 
                           avail * p_t_tilde2 * (1-p_t_tilde2) * (exp_Xdm_beta2^(-1) * as.vector(EY2) - as.vector(EY0)))* Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = 2 * p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside estimator 2():")
      message(cond)
      return(list(root = rep(NaN, 2 * p), msg = cond,
                  f.root = rep(NaN, 2 * p)))
    })
  
  beta_hat <- matrix(solution$root, nrow = p)
  
  ### asymptotic variance ###
  
  ### Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
  
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, 2 * p, 2 * p))
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, 2 * p, 2 * p))

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(t(Xdm[it, ]) %*% beta_hat)
    }

    
    pre_multiplier <- exp(- A1[it] * Xbeta[1] - A2[it] * Xbeta[2]) * WCLS_weight[it]
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow =  2 * p, ncol = 2 * p)
    partialD_partialtheta[1:p, 1:p] <- -pre_multiplier * A1[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[1:p, (p + 1):(2*p)] <- -pre_multiplier * A2[it] * cA_tilde1[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(p + 1):(2*p), 1:p] <- -pre_multiplier * A1[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(p + 1):(2*p), (p + 1):(2*p)] <- -pre_multiplier * A2[it] * cA_tilde2[it] * (Xdm[it, ] %o% Xdm[it, ])

    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - A1[it] - A2[it]) + EY1[it] * A1[it] + EY2[it] * A2[it]
    r_term <- (Y[it] - EH) * avail[it]

    # # D_term = D^{(t),T} (dim = (p) * 1)
    D_term <- pre_multiplier * c(cA_tilde1[it] * Xdm[it, ], cA_tilde2[it] * Xdm[it, ])

    # # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- rep(0, 2*p)
    # partialr_partialtheta_collected[it, ] <- partialr_partialtheta

    # ####add double robust term### p * 1
    add_term <- avail[it] * c(p_t_tilde1[it] * (1 - p_t_tilde1[it]) * (EY1[it] * exp(-Xbeta[1]) - EY0[it]) * Xdm[it, ], p_t_tilde2[it] * (1 - p_t_tilde2[it]) * (EY2[it] * exp(-Xbeta[2]) - EY0[it]) * Xdm[it, ])
    add_term_partial <- matrix(0, nrow =  2 * p, ncol = 2 * p)
    add_term_partial[1:p, 1:p] <- -avail[it] * p_t_tilde1[it] * (1 - p_t_tilde1[it]) * EY1[it] * exp(-Xbeta[1]) * (Xdm[it, ] %o% Xdm[it, ])
    add_term_partial[(p + 1):(2*p), (p + 1):(2*p)] <- -avail[it] * p_t_tilde2[it] * (1 - p_t_tilde2[it]) * EY2[it] * exp(-Xbeta[2]) * (Xdm[it, ] %o% Xdm[it, ])

    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta + add_term_partial
    Sigman_summand[it, , ] <-  (D_term * r_term + add_term) %*% (r_term * t(D_term) + t(add_term))
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size / time_length
  Mn_inv <- solve(Mn)
  
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size / time_length
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size / time_length
  beta_se <- sqrt(diag(varcov)[1:(2 * p)])
  
  
  beta_hat = as.vector(beta_hat)
  names(beta_hat) <- names(beta_se) <- rep(Xnames, 2)
  
  return(list(beta_hat = beta_hat, beta_se = beta_se,
              varcov = varcov,
              dims = list(p = p),
              f.root = solution$f.root))
}
