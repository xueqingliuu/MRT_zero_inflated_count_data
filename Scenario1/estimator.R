library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package

get_alpha_beta_from_multiroot_result <- function(root, p, q)
{
  if (p == 1) {
    beta_root <- root$root[q+1]
  } else {
    beta_root <- as.matrix(root$root[(q+1) : (q+p)])
  }
  if (q == 1) {
    alpha_root <- root$root[1]
  } else {
    alpha_root <- as.matrix(root$root[1:q])
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
    
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 1) {
      p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
      p_t_tilde <- rand_prob_tilde
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde <- A - p_t_tilde
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Intercept", control_varname)

  
  WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  

  ### estimation ###
  
  estimating_equation <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])
    
    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1)
    
    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * WCLS_weight * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde * Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p + q)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, p + q), msg = cond,
                  f.root = rep(NaN, p + q)))
    })
  
  estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
  alpha_hat <- as.vector(estimator$alpha)
  beta_hat <- as.vector(estimator$beta)
  
  ### 3. asymptotic variance ###
  
  ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
  
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
  r_term_collected <- rep(NA, total_person_decisionpoint)
  D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_hat
    } else {
      Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
    }
    
    pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight[it]
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
    partialD_partialtheta[1:q, 1:q] <- 0
    partialD_partialtheta[1:q, (q+1):(q+p)] <- - pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
    partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
    r_term_collected[it] <- r_term
    
    # D_term = D^{(t),T} (dim = (p+q) * 1)
    D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde[it] * Xdm[it, ])
    D_term_collected[, it] <- D_term
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
 
  Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    
    Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  alpha_se <- sqrt(diag(varcov)[1:q])
  beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
  

  names(alpha_hat) <- names(alpha_se) <- Znames
  names(beta_hat) <- names(beta_se) <- Xnames
  
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
  
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 1) {
      p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
      p_t_tilde <- rand_prob_tilde
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  ###estimation: nuisance parameter  ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  
  
  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  
  ### estimation: parameter of interest [beta] ###
  
  estimating_equation <- function(beta) {
    
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - p_t_tilde) + as.vector(EY1) * p_t_tilde * exp_Xdm_beta^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1)
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde * Xdm[, i] )
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- solution$root
  
  ### asymptotic variance ###
  
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }

    pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight[it]
    
    partialD_partialtheta <- matrix(NA, nrow = p, ncol = p)
    partialD_partialtheta[1:p, 1:p] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - p_t_tilde[it]) + EY1[it] * p_t_tilde[it] * exp(-Xbeta)
    r_term <- (Y[it] - EH * exp(A[it] * Xbeta)) * avail[it]
    
    # D_term = D^{(t),T} (dim = (p) * 1)
    D_term <- pre_multiplier * cA_tilde[it] * Xdm[it, ]

    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- EY1[it] * p_t_tilde[it] * exp(-Xbeta) * exp(A[it] * Xbeta) * avail[it] * Xdm[it, ]  - EH * exp(A[it] * Xbeta) * A[it] * Xdm[it, ] * avail[it]

    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    Sigman_summand[it, , ] <- (D_term * r_term) %*% (r_term * t(D_term))
  }
  
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size/ time_length
  Mn_inv <- solve(Mn)
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size/ time_length
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size/ time_length
  beta_se <- sqrt(diag(varcov)[1:p])

  
  names(beta_hat) <- names(beta_se) <- Xnames
  
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
   
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 1) {
      p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
      p_t_tilde <- rand_prob_tilde
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))
  
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  ### estimation: nuisance parameter [] ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  
  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  
  ### estimation: parameter of interest [beta] ###
  
  estimating_equation <- function(beta) {
    
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - A) + as.vector(EY1) * A * exp_Xdm_beta^(-1)
    # exp_Zdm_alpha <- as.vector(EY0) * (1 - p_t_tilde) + as.vector(EY1) * p_t_tilde * exp_Xdm_beta^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_AXdm_beta^(-1) 

    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( (weight * residual * avail * WCLS_weight * cA_tilde + 
                       avail * p_t_tilde * (1-p_t_tilde) * (exp_Xdm_beta^(-1) * as.vector(EY1) - as.vector(EY0)))* Xdm[, i] )
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- solution$root
  
  ### asymptotic variance ###
  
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }

    
    pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight[it]
    
    partialD_partialtheta <- matrix(NA, nrow = p, ncol = p)
    partialD_partialtheta[1:p, 1:p] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - A[it]) + EY1[it] * A[it]
    r_term <- (Y[it] - EH ) * avail[it]

    D_term <- pre_multiplier * cA_tilde[it] * Xdm[it, ]

    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <-  rep(0, p)
    
    
    ####add double robust term### p * 1
    add_term <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * (EY1[it] * exp(-Xbeta) - EY0[it]) * Xdm[it, ]
    add_term_partial <- -avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * EY1[it] * exp(-Xbeta) * (Xdm[it, ] %o% Xdm[it, ])
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta + add_term_partial
    Sigman_summand[it, , ] <-  (D_term * r_term + add_term) %*% (r_term * t(D_term) + t(add_term))
  }
  
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size / time_length
  Mn_inv <- solve(Mn)
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size / time_length
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size / time_length
  beta_se <- sqrt(diag(varcov)[1:p])

  

  names(beta_hat) <- names(beta_se) <- Xnames
  
  return(list(beta_hat = beta_hat, beta_se = beta_se, 
              varcov = varcov,
              dims = list(p = p),
              f.root = solution$f.root))
}

log_linear_GEE_geepack <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    estimator_initial_value = NULL,
    corstr = "independence" 
){
  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Intercept", control_varname)
  
  control_summed <- paste0(control_varname, collapse = " + ")
  if (p > 1) {
    moderator_summed <- paste0("* (", paste0(moderator_varname, collapse = " + "), ")")
  } else {
    moderator_summed <- ""
  }
  
  gee_formula <- as.formula(paste(outcome_varname, "~", control_summed, "+", treatment_varname, moderator_summed))
  fit_geepack <- geeglm(gee_formula, data = dta, corstr = corstr, id = dta[, id_varname], family = poisson("log"))
  
  
  alpha_hat <- fit_geepack$coefficients[1:q]
  beta_hat <- fit_geepack$coefficients[(q+1):(q+p)]
  varcov <- vcov(fit_geepack)
  alpha_se <- sqrt(diag(varcov)[1:q])
  beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])

  
  names(alpha_hat) <- names(alpha_se) <- Znames
  names(beta_hat) <- names(beta_se) <- Xnames
  
  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = rep(-1, p+q)))
}


fit_ECE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
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
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  

  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Intercept", control_varname)
  
  ### 2. estimation ###
  
  estimating_equation <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])
    
    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_negAXdm_beta <- exp_AXdm_beta^(-1)
    exp_negXdm_beta <- exp_Xdm_beta^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- - exp_negAXdm_beta / ( p_t + exp_negXdm_beta * (1 - p_t) )
    
    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * cA * Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p + q)
  }
  
  # browser()
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee():")
      message(cond)
      return(list(root = rep(NaN, p + q), msg = cond,
                  f.root = rep(NaN, p + q)))
    })
  
  estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
  alpha_hat <- as.vector(estimator$alpha)
  beta_hat <- as.vector(estimator$beta)
  
  ###  asymptotic variance ###

  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
  r_term_collected <- rep(NA, total_person_decisionpoint)
  D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_hat
    } else {
      Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
    }
  
    denom <- p_t[it] + exp(-Xbeta) * (1 - p_t[it])
    W1 <- exp(- A[it] * Xbeta) / (denom^2)
    W3 <- A[it] * denom - exp(-Xbeta) * (1 - p_t[it])  
    
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
    partialD_partialtheta[1:q, 1:q] <- 0
    partialD_partialtheta[1:q, (q+1):(q+p)] <- W1 * W3 * (Zdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
    partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
    r_term_collected[it] <- r_term
    
    # D_term = D^{(t),T} (vector of length (p+q))
    D_term <- -exp(- A[it] * Xbeta) / ( p_t[it] + exp(-Xbeta)  * (1 - p_t[it]) ) *
      c(Zdm[it, ], cA[it] * Xdm[it, ])
    D_term_collected[, it] <- D_term
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
    partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  

  Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  alpha_se <- sqrt(diag(varcov)[1:q])
  beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
  
  
  names(alpha_hat) <- names(alpha_se) <- Znames
  names(beta_hat) <- names(beta_se) <- Xnames
  
  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}

fit_ECE_NonP <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
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
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  

  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  
  ### estimation: nuisance parameter ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  
  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  ### estimation equation###
  
  estimating_equation <- function(beta) {
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - p_t) + as.vector(EY1) * p_t * exp_Xdm_beta^(-1)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_negAXdm_beta <- exp_AXdm_beta^(-1)
    exp_negXdm_beta <- exp_Xdm_beta^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- - exp_negAXdm_beta / ( p_t + exp_negXdm_beta * (1 - p_t) )
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( weight * residual * avail * cA * Xdm[, i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- solution$root
  
  ### asymptotic variance ###
  
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }

    denom <- p_t[it] + exp(-Xbeta) * (1 - p_t[it])
    W1 <- exp(- A[it] * Xbeta) / (denom^2)
    W3 <- A[it] * denom - exp(-Xbeta) * (1 - p_t[it])
    
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow = p, ncol = p)
    partialD_partialtheta[(1):(p), (1):(p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - p_t[it]) + EY1[it] * p_t[it] * exp(-Xbeta)
    r_term <- (Y[it] - EH * exp(A[it] * Xbeta)) * avail[it]
    
    # D_term = D^{(t),T} (vector of length (p+q))
    D_term <- - exp(- A[it] * Xbeta) / ( p_t[it] + exp(-Xbeta) * (1 - p_t[it]) ) * cA[it] * Xdm[it, ]
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
    partialr_partialtheta <- - EH * exp(A[it] * Xbeta) * A[it] * Xdm[it, ] * avail[it]
    
    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    Sigman_summand[it, , ] <- (D_term * r_term) %*% (r_term * t(D_term))
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size / time_length
  Mn_inv <- solve(Mn)
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size / time_length
   
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size / time_length
  beta_se <- sqrt(diag(varcov)[(1):(p)])

  
  names(beta_hat) <- names(beta_se) <- Xnames
  
  return(list(beta_hat = beta_hat, 
              beta_se = beta_se, 
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}

fit_G_est_Yu <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
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
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
  

  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Y","Intercept", control_varname)
  
  
  ### estimation: nuisance parameter ###
  df4 = data.frame(Y=as.vector(Y), Zdm)
  colnames(df4) = Znames
  df40 = df4[A==0,]
  df40_prob = df40
  df40_prob$ind = df40$Y > 0
  df40_mean = df40[which(df40_prob$Y > 0),]
  df41 = df4[A==1,]
  df41_prob = df41
  df41_prob$ind = df41$Y > 0
  df41_mean = df41[which(df41_prob$Y > 0),]
  
  formula1 <- as.formula(paste0("ind", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  formula2 <- as.formula(paste0("Y", "~", "s(",paste0(control_varname, collapse = "+"), ",k=3)"))
  fit40_prob = gam(formula1, data=df40_prob, family = poisson(link='log'))
  fit40_mean = gam(formula2, data=df40_mean, family = gaussian(link="log"))
  fit41_prob = gam(formula1, data=df41_prob, family = poisson(link='log'))
  fit41_mean = gam(formula2, data=df41_mean, family = gaussian(link="log"))
  
  EY1 = predict(fit41_prob, newdata=df4, type = "response") * predict(fit41_mean, newdata=df4, type = "response")
  EY0 = predict(fit40_prob, newdata=df4, type = "response") * predict(fit40_mean, newdata=df4, type = "response")
  
  ### estimation equation###
  
  estimating_equation <- function(beta) {
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_Zdm_alpha <- as.vector(EY0) * (1 - p_t) + as.vector(EY1) * p_t * exp_Xdm_beta^(-1)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_negAXdm_beta <- exp_AXdm_beta^(-1)
    exp_negXdm_beta <- exp_Xdm_beta^(-1)
    
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- exp_negAXdm_beta
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( weight * residual * avail * cA * Xdm[, i])
    }
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- solution$root
  
  ### asymptotic variance ###
  
  Sigman_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  
  for (it in 1:total_person_decisionpoint) {  
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.  
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    
    W1 <- exp(- A[it] * Xbeta)
    W3 <- - A[it]
    
    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow = p, ncol = p)
    partialD_partialtheta[(1):(p), (1):(p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
    
    # r_term = r^(t) (scalar)
    EH <- EY0[it] * (1 - p_t[it]) + EY1[it] * p_t[it] * exp(-Xbeta)
    r_term <- (Y[it] - EH * exp(A[it] * Xbeta)) * avail[it]

    # D_term = D^{(t),T} (vector of length (p+q))
    D_term <- exp(- A[it] * Xbeta) * cA[it] * Xdm[it, ]

    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
    partialr_partialtheta <- - EH * exp(A[it] * Xbeta) * A[it] * Xdm[it, ] * avail[it]

    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    Sigman_summand[it, , ] <- (D_term * r_term) %*% (r_term * t(D_term))
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size / time_length
  Mn_inv <- solve(Mn)
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size / time_length
   
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size / time_length
  beta_se <- sqrt(diag(varcov)[(1):(p)])

  names(beta_hat) <- names(beta_se) <- Xnames
  
  return(list(beta_hat = beta_hat, 
              beta_se = beta_se, 
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}
