
expit <- function(x){
    return(exp(x)/(1+exp(x)))
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
        prob_a <- expit(eta1*A + eta2*dta$S[row_index])
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

