rm(list = ls())

source("estimator_observational.R")
source("data_generation.R")
source("other_functions.R")

library(tidyverse)
library(MASS)
library(foreach)
library(doMC)
library(doRNG)
library(mgcv)
max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))

sample_size <- 100
total_Ts <- c(30, 100, 150)
nsim <- 1000

set.seed(1)
control_vars <- "S"
moderator_vars <- "S"

result_df_collected <- data.frame()

for (i_ss in 1:length(total_Ts)) {
  
  total_T <- total_Ts[i_ss]
  
  writeLines(c(""), "~/Downloads/log.txt")
  sink("~/Downloads/log.txt", append=FALSE)
  result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
    if (isim %% 10 == 0) {
      cat(paste("Starting iteration",isim,"\n"))
    }
      dta <- generate_data(sample_size, total_T)
      rand_prob_tilde = mean(dta$A)
      
      EMEE <- fit_EMEE(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        rand_prob = rand_prob_tilde,
        rand_prob_tilde_varname = NULL,
        rand_prob_tilde = rand_prob_tilde,
        estimator_initial_value = NULL
      )
      
      EMEE_NonP <- fit_EMEE_NonP(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        rand_prob = rand_prob_tilde,
        rand_prob_tilde_varname = NULL, 
        rand_prob_tilde = rand_prob_tilde,        
        estimator_initial_value = NULL
      )
      
      
      ECE <- fit_ECE(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        rand_prob = rand_prob_tilde,
        estimator_initial_value = NULL
      )
      
      ECE_NonP <- fit_ECE_NonP(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        rand_prob = rand_prob_tilde,
        estimator_initial_value = NULL
      )
      
      DR_EMEE_NonP <- fit_DR_EMEE_NonP(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        rand_prob = rand_prob_tilde,
        rand_prob_tilde_varname = NULL, 
        rand_prob_tilde = rand_prob_tilde,        
        estimator_initial_value = NULL
      )
      
      GEE_IND <- log_linear_GEE_geepack(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        estimator_initial_value = NULL,
        corstr = "independence"
      )
      
      GEE_EXCH <- log_linear_GEE_geepack(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = control_vars,
        moderator_varname = moderator_vars,
        estimator_initial_value = NULL,
        corstr = "exchangeable"
      )
      output <- list(list(ECE = ECE, ECE_NonP = ECE_NonP, EMEE = EMEE, 
                          EMEE_NonP = EMEE_NonP, DR_EMEE_NonP = DR_EMEE_NonP, 
                          GEE_IND = GEE_IND, GEE_EXCH = GEE_EXCH))
  }
  
  
  ee_names <- c("ECE", "ECE_NonP", "EMEE", "EMEE_NonP", "DR_EMEE_NonP", "GEE_IND", "GEE_EXCH")
  beta_names <- c("Intercept", moderator_vars)
  num_estimator <- length(ee_names)
  
  

  beta <- simplify2array(lapply(result, function(l) matrix(c(l$ECE$beta_hat, l$ECE_NonP$beta_hat, l$EMEE$beta_hat, l$EMEE_NonP$beta_hat,
                                                             l$DR_EMEE_NonP$beta_hat, l$GEE_IND$beta_hat, l$GEE_EXCH$beta_hat),
                                                           nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
  beta_se <- simplify2array(lapply(result, function(l) matrix(c(l$ECE$beta_se, l$ECE_NonP$beta_se, l$EMEE$beta_se, l$EMEE_NonP$beta_se, 
                                                                l$DR_EMEE_NonP$beta_se, l$GEE_IND$beta_se, l$GEE_EXCH$beta_se),
                                                              nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
  
  result <- compute_result_beta(beta_true, beta, beta_se, moderator_vars, control_vars, significance_level = 0.05)
  result_df <- data.frame(time = rep(total_T, num_estimator),
                          est = ee_names,
                          bias1 = result$bias[,1],
                          sd1 = result$sd[,1],
                          rmse1 = result$rmse[,1],
                          cp1 = result$coverage_prob[,1],
                          se1 = result$se[,1],
                          bias2 = result$bias[,2],
                          sd2 = result$sd[,2],
                          rmse2 = result$rmse[,2],
                          cp2 = result$coverage_prob[,2],
                          se2 = result$se[,2])
  names(result_df) <- c("ss", "est", "bias1", "sd1", "rmse1", "cp1", "se1",
                        "bias2", "sd2", "rmse2", "cp2", "se2")
  rownames(result_df) <- NULL
  
  result_df_collected <- rbind(result_df_collected, result_df)
}

saveRDS(result_df_collected, file = "result_simulation_mod_zinb.RDS")

