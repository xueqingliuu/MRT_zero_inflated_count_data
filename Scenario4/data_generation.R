
expit <- function(x){
    return(exp(x)/(1+exp(x)))
}

generate_data <- function(sample_size, total_T) {
    
    baseline_Y_S0 <- 2.2
    baseline_Y_S1 <- 2.5
    baseline_Y_S2 <- 2.4
    
    beta_10 <- 0.1
    beta_11 <- 0.3
    beta_20 <- 0.1
    beta_21 <- 0.1
    
    eta1 <- -0.2
    eta2 <- -0.3
    eta3 <- 0.5
    
    gam1 <- -0.4
    gam2 <- 0.1
    gam3 <- 0.1
    
    
    df_names <- c("userid", "day", "Y", "A1","A2", "S", "S2", "prob_zero","prob_Y", "prob_Y_A0", "prob_A1", "prob_A2")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    A1 <- rep(0, sample_size)
    A2 <- rep(0, sample_size)
    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
    
        dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
        dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0)
        prob_a <- expit(eta1*A1 + eta2*A2 + eta3*dta$S[row_index])
        prob_a1 <- prob_a * 1/2
        prob_a2 <- prob_a * 1/2
        
        # prob_a <- rep(0.2, sample_size)
        dta$prob_A1[row_index] <- prob_a1
        dta$prob_A2[row_index] <- prob_a2
        A <- c()
        for (i in 1:sample_size) {
          A[i] <- sample(c(0,1,2), 1, replace = T, prob = c(1-prob_a1[i]-prob_a2[i], prob_a1[i], prob_a2[i]))
        }
        

       dta$A1[row_index] <- ifelse(A==1, 1, 0)
       dta$A2[row_index] <- ifelse(A==2, 1, 0)
        

        ### ZINB ###
        dta$prob_Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                           ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
        dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A1[row_index] * (beta_10 + beta_11 * dta$S[row_index]) + dta$A2[row_index] * (beta_20 + beta_21 * dta$S[row_index]))
        dta$prob_zero[row_index] <- exp(gam1* (dta$S[row_index]+0.1) + gam2 * dta$S[row_index] * dta$A1[row_index] + gam3 * dta$S[row_index] * dta$A2[row_index])
        
        zero <- rbinom(sample_size, 1, dta$prob_zero[row_index])
        dta$Y[row_index] <- zero * rnbinom(sample_size, size = 1, mu = dta$prob_Y[row_index])
        
    }
    
    return(dta)
}

## true beta for (Intercept, S)
beta_true <- c(0.1, 0.4, 0.1, 0.2)

## true beta for Intercept
beta_true_marginal <- c(0.4599, 0.2672)



# compute marginal beta_true
if (0) {
    
    ### analytically, we have ###
    
    baseline_Y_S0 <- 2.2
    baseline_Y_S1 <- 2.5
    baseline_Y_S2 <- 2.4
    
    beta_10 <- 0.1
    beta_11 <- 0.3
    beta_20 <- 0.1
    beta_21 <- 0.1
    
    # prob_a <- 0.2
    
    gam1 <- -0.4
    gam2 <- 0.1
    gam3 <- 0.1

    numerator1 <- exp(gam1 * 0.1) * baseline_Y_S0 * exp(beta_10) + exp(gam1 * 1.1 + gam2) * baseline_Y_S1 * exp(beta_10 + beta_11) + exp(gam1 * 2.1 + 2 * gam2) * baseline_Y_S2 * exp(beta_10 + 2 * beta_11)
    denominator <- exp(gam1 * 0.1) * baseline_Y_S0 + exp(gam1 * 1.1) * baseline_Y_S1 + exp(gam1 * 2.1) * baseline_Y_S2
    log(numerator1 / denominator)
    
    numerator2 <- exp(gam1 * 0.1) * baseline_Y_S0 * exp(beta_20) + exp(gam1 * 1.1 + gam2) * baseline_Y_S1 * exp(beta_20 + beta_21) + exp(gam1 * 2.1 + 2 * gam2) * baseline_Y_S2 * exp(beta_20 + 2 * beta_21)
    log(numerator2 / denominator)

}



