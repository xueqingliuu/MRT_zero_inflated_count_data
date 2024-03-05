#data analysis of drink less MRT (count outcome)
library(MASS)
library(parallel)
source("estimator.R")
TQ_session_modules <- readRDS("TQ_session_modules.rds")
# TQ_session_length <- readRDS("TQ_session_length.rds")
a <- readRDS("FINAL Dataset_A.rds")

num_users <- length(unique(TQ_session_modules$ID)) #349
num_days <- length(unique(TQ_session_modules$days_since_download)) #31 ?there is 0?

TQ_session_modules = TQ_session_modules[TQ_session_modules$days_since_download!=0,]
num_days <- length(unique(TQ_session_modules$days_since_download)) #30

id <- unique(TQ_session_modules$ID)
ID <- rep(unique(TQ_session_modules$ID), each = num_days)
Day <- rep(c(1:num_days), num_users)

number_of_screen1 <- c() 
yesterday_screen_view <- c() 
find_screen_number1 <- function(temp2) {
  num <- nrow(temp2)
  temp2$screen_view <- format(as.POSIXct(temp2$screen_view), format = "%H:%M:%S")
  # if (num == 1) {
  #   count = ifelse(is.na(temp2$screen_view), 0, 1)
  # } else {
  count = NROW(temp2[!is.na(temp2$screen_view)&temp2$screen_view >= "20:00:00" & temp2$screen_view <= "21:00:00",])
  # }
  return(count)
}

for (i in (1:num_users)) {
  temp1 <- TQ_session_modules[TQ_session_modules$ID==id[i],]
  for (j in (1:num_days)) {
    temp2 <- temp1[temp1$days_since_download==j,]
    num_views <- find_screen_number1(temp2)
    number_of_screen1[(i-1)*num_days+j] <- num_views
    if (j < num_days) {
      yesterday_screen_view[(i-1)*num_days + j + 1] <- num_views
    }
    
  }
}


#intervention
library(dplyr)
TQ_session_modules <- TQ_session_modules %>% distinct(ID, days_since_download, .keep_all = TRUE)
Intervention <- TQ_session_modules$message
levels(Intervention)[1] <- "No notification"
levels(Intervention)[18] <- "Standard notification"
levels(Intervention)[-c(1,18)] <- "New notification"
#encoding with dummy variables

#covariates
#employment type: 0 non-manual, 1 manual, 2 other
#gender: 0 male 1 female
dataset <- data.frame(ID, Day, number_of_screen1, yesterday_screen_view, Intervention, 
                      "gender" = TQ_session_modules$gender, "days_since_download" = TQ_session_modules$days_since_download, 
                      "AUDIT_score" = TQ_session_modules$AUDIT_score, "age" = TQ_session_modules$age,
                      "employment_type" = TQ_session_modules$employment_type)

library(fastDummies)
dataset <- dummy_cols(dataset, select_columns = c( "Intervention", "employment_type"), remove_first_dummy = T)
colnames(dataset)[13] <- "Manual_employ"
colnames(dataset)[14] <- "Other_employ"
colnames(dataset)[11] <- "New_notification"
colnames(dataset)[12] <- "Standard_notification"
dataset$age <- as.numeric(scale(dataset$age))
dataset$treatment <- ifelse(Intervention=="No notification", 0, 1)

dataset$AUDIT_score <- as.numeric(scale(dataset$AUDIT_score))

dataset$yesterday_screen_view[is.na(dataset$yesterday_screen_view)] <- 0

prob_A <- ifelse(dataset$treatment==1, 0.6, 0.4)
dataset$prob_A <- prob_A


library(mgcv)

######Primary: give vs not give, S_t = 1#######


control_vars = c('age','AUDIT_score','Day','gender','Manual_employ', 'Other_employ', 'yesterday_screen_view')
moderator_vars = NULL

fit_wcls <- weighted_centered_least_square(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.6,
  estimator_initial_value = NULL
)

fit_wcls$beta_hat # 1.107017 
fit_wcls$beta_se # 0.1259011

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_wcls <- round(c(rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se)),3)
# (0.8593891,  1.3546449)
p_wcls <- 2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_wcls_use <- paste("(", paste(ci_wcls[1],ci_wcls[2], sep = ","), ")", sep = "")


fit_nonp_wcls <- two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)

fit_nonp_wcls$beta_hat # 1.119966  
fit_nonp_wcls$beta_se # 0.09425942 

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_nonp_wcls <- round(c(rbind(fit_nonp_wcls$beta_hat - t_quantile * fit_nonp_wcls$beta_se,
      fit_nonp_wcls$beta_hat + t_quantile * fit_nonp_wcls$beta_se)), 3)
# (0.9273521,  1.3002999)
p_nonp_wcls <- 2 * pt(abs(fit_nonp_wcls$beta_hat) / fit_nonp_wcls$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_nonp_wcls_use <- paste("(", paste(ci_nonp_wcls[1],ci_nonp_wcls[2], sep = ","), ")", sep = "")



fit_dr_nonp_wcls <- DR_two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)

fit_dr_nonp_wcls$beta_hat # 1.119966  
fit_dr_nonp_wcls$beta_se # 0.09425942 

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_dr_nonp_wcls <- round(c(rbind(fit_dr_nonp_wcls$beta_hat - t_quantile * fit_dr_nonp_wcls$beta_se,
                                 fit_dr_nonp_wcls$beta_hat + t_quantile * fit_dr_nonp_wcls$beta_se)), 3)
# (0.9273521,  1.3002999)
p_nonp_wcls <- 2 * pt(abs(fit_dr_nonp_wcls$beta_hat) / fit_dr_nonp_wcls$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_dr_nonp_wcls_use <- paste("(", paste(ci_dr_nonp_wcls[1],ci_dr_nonp_wcls[2], sep = ","), ")", sep = "")



##### Primary: three treatments, S_t=1#####
source("estimator_three_cate.R")
prob_A1 <- prob_A2 <- rep(0.3, nrow(dataset))
dataset$prob_A1 <- prob_A1
dataset$prob_A2 <- prob_A2

fit_wcls <- weighted_centered_least_square(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = c("Standard_notification", "New_notification"),
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = c("prob_A1", "prob_A2"),
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = c(0.3, 0.3),
  estimator_initial_value = NULL
)

fit_wcls$beta_hat # 1.241604   0.965
fit_wcls$beta_se # 0.1259011

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_wcls <- round(c(rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se,
                         fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se)),3)
ci# (0.8593891,  1.3546449)
p_wcls <- 2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_wcls_use1 <- paste("(", paste(ci_wcls[1], ci_wcls[2], sep = ","), ")", sep = "")
ci_wcls_use2 <- paste("(", paste(ci_wcls[3], ci_wcls[4], sep = ","), ")", sep = "")



fit_nonp_wcls <- two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = c("Standard_notification", "New_notification"),
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = c("prob_A1", "prob_A2"),
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = c(0.3, 0.3),
  estimator_initial_value = NULL
)
fit_nonp_wcls$beta_hat # 1.119966  
fit_nonp_wcls$beta_se # 0.09425942 

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_nonp_wcls <- round(c(rbind(fit_nonp_wcls$beta_hat - t_quantile * fit_nonp_wcls$beta_se,
                              fit_nonp_wcls$beta_hat + t_quantile * fit_nonp_wcls$beta_se)), 3)
# (0.9273521,  1.3002999)
p_nonp_wcls <- 2 * pt(abs(fit_nonp_wcls$beta_hat) / fit_nonp_wcls$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27


ci_nonp_wcls_use1 <- paste("(", paste(ci_nonp_wcls[1], ci_nonp_wcls[2], sep = ","), ")", sep = "")
ci_nonp_wcls_use2 <- paste("(", paste(ci_nonp_wcls[3], ci_nonp_wcls[4], sep = ","), ")", sep = "")



fit_dr_nonp_wcls <- DR_two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = c("Standard_notification", "New_notification"),
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = c("prob_A1", "prob_A2"),
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = c(0.3, 0.3),
  estimator_initial_value = NULL
)
fit_dr_nonp_wcls$beta_hat # 1.119966  
fit_dr_nonp_wcls$beta_se # 0.09425942 

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_dr_nonp_wcls <- round(c(rbind(fit_dr_nonp_wcls$beta_hat - t_quantile * fit_dr_nonp_wcls$beta_se,
                              fit_dr_nonp_wcls$beta_hat + t_quantile * fit_dr_nonp_wcls$beta_se)), 3)
# (0.9273521,  1.3002999)
p_dr_nonp_wcls <- 2 * pnorm((fit_dr_nonp_wcls$beta_hat) / fit_dr_nonp_wcls$beta_se,  lower.tail = FALSE)
# 4.791444e-27


ci_nonp_wcls_use1 <- paste("(", paste(ci_nonp_wcls[1], ci_nonp_wcls[2], sep = ","), ")", sep = "")
ci_nonp_wcls_use2 <- paste("(", paste(ci_nonp_wcls[3], ci_nonp_wcls[4], sep = ","), ")", sep = "")


##### Secondary:  S_t = 1 + Day#####
source("estimator.R")
control_vars = c('age','AUDIT_score','Day','gender','Manual_employ', 'Other_employ', 'yesterday_screen_view')
moderator_vars = 'Day'

fit_wclsmod <- weighted_centered_least_square(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.6,
  estimator_initial_value = NULL
)

fit_wclsmod$beta_hat # 1.32914887    -0.01719877
fit_wclsmod$beta_se # 0.17711650    0.01067966

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_wclsmod <- round(c(rbind(fit_wclsmod$beta_hat - t_quantile * fit_wclsmod$beta_se,
      fit_wclsmod$beta_hat + t_quantile * fit_wclsmod$beta_se)),3)
# (0.9807884,  1.6775094)    (-0.03820399, 0.00380646)
p_wclsmod <- 2 * pnorm(abs(fit_wclsmod$beta_hat ) / fit_wclsmod$beta_se_adjusted, lower.tail = FALSE)
# 5.269673e-13    1.082175e-01

ci_wclsmod_use1 <- paste("(", paste(ci_wclsmod[1], ci_wclsmod[2], sep = ","), ")", sep = "")
ci_wclsmod_use2 <- paste("(", paste(ci_wclsmod[3], ci_wclsmod[4], sep = ","), ")", sep = "")

fit_dr_nonp_wclsmod <- DR_two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)
fit_dr_nonp_wclsmod$beta_hat # 1.49933525   -0.02946211
fit_dr_nonp_wclsmod$beta_se # 0.153213587  0.008847506

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)  
ci_dr_nonp_wclsmod <- round(c(rbind(fit_dr_nonp_wclsmod$beta_hat - t_quantile * fit_dr_nonp_wclsmod$beta_se,
                                 fit_dr_nonp_wclsmod$beta_hat + t_quantile * fit_dr_nonp_wclsmod$beta_se)),3)
# (1.197988,  1.800682)   (-0.04686377, -0.01206044)
p_dr_nonp_wclsmod <- 2 * pnorm(abs(fit_dr_nonp_wclsmod$beta_hat) / fit_dr_nonp_wclsmod$beta_se, lower.tail = FALSE)
# 3.998510e-20   9.620942e-04

ci_dr_nonp_wclsmod_use1 <- paste("(", paste(ci_dr_nonp_wclsmod[1], ci_dr_nonp_wclsmod[2], sep = ","), ")", sep = "")
ci_dr_nonp_wclsmod_use2 <- paste("(", paste(ci_dr_nonp_wclsmod[3], ci_dr_nonp_wclsmod[4], sep = ","), ")", sep = "")


fit_nonp_wclsmod <- two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)
fit_nonp_wclsmod$beta_hat # 1.49933525   -0.02946211
fit_nonp_wclsmod$beta_se # 0.153213587  0.008847506

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)  
ci_nonp_wclsmod <- round(c(rbind(fit_nonp_wclsmod$beta_hat - t_quantile * fit_nonp_wclsmod$beta_se,
                                 fit_nonp_wclsmod$beta_hat + t_quantile * fit_nonp_wclsmod$beta_se)),3)
# (1.197988,  1.800682)   (-0.04686377, -0.01206044)
p_nonp_wclsmod <- 2 * pnorm(abs(fit_nonp_wclsmod$beta_hat) / fit_nonp_wclsmod$beta_se,  lower.tail = FALSE)
# 3.998510e-20   9.620942e-04

ci_nonp_wclsmod_use1 <- paste("(", paste(ci_nonp_wclsmod[1], ci_nonp_wclsmod[2], sep = ","), ")", sep = "")
ci_nonp_wclsmod_use2 <- paste("(", paste(ci_nonp_wclsmod[3], ci_nonp_wclsmod[4], sep = ","), ")", sep = "")



#####S_t = 1 + Yesterday#####
control_vars = c('age','AUDIT_score','Day','gender','Manual_employ', 'Other_employ', 'yesterday_screen_view')
moderator_vars = 'yesterday_screen_view'

fit_wclsmod <- weighted_centered_least_square(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.6,
  estimator_initial_value = NULL
)

fit_wclsmod$beta_hat # 1.32914887    -0.01719877
fit_wclsmod$beta_se # 0.17711650    0.01067966

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)  
ci_wclsmod <- round(c(rbind(fit_wclsmod$beta_hat - t_quantile * fit_wclsmod$beta_se,
                            fit_wclsmod$beta_hat + t_quantile * fit_wclsmod$beta_se)),3)
# (0.9807884,  1.6775094)    (-0.03820399, 0.00380646)
p_wclsmod <- 2 * pnorm(abs(fit_wclsmod$beta_hat) / fit_wclsmod$beta_se, lower.tail = FALSE)
# 5.269673e-13    1.082175e-01

ci_wclsmod_use1 <- paste("(", paste(ci_wclsmod[1], ci_wclsmod[2], sep = ","), ")", sep = "")
ci_wclsmod_use2 <- paste("(", paste(ci_wclsmod[3], ci_wclsmod[4], sep = ","), ")", sep = "")

fit_nonp_wclsmod <- two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)
fit_nonp_wclsmod$beta_hat # 1.49933525   -0.02946211
fit_nonp_wclsmod$beta_se # 0.153213587  0.008847506

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)   
ci_nonp_wclsmod <- round(c(rbind(fit_nonp_wclsmod$beta_hat - t_quantile * fit_nonp_wclsmod$beta_se,
                                 fit_nonp_wclsmod$beta_hat + t_quantile * fit_nonp_wclsmod$beta_se)),3)
# (1.197988,  1.800682)   (-0.04686377, -0.01206044)
p_nonp_wclsmod <- 2 * pnorm(abs(fit_nonp_wclsmod$beta_hat) / fit_nonp_wclsmod$beta_se,  lower.tail = FALSE)
# 3.998510e-20   9.620942e-04

ci_nonp_wclsmod_use1 <- paste("(", paste(ci_nonp_wclsmod[1], ci_nonp_wclsmod[2], sep = ","), ")", sep = "")
ci_nonp_wclsmod_use2 <- paste("(", paste(ci_nonp_wclsmod[3], ci_nonp_wclsmod[4], sep = ","), ")", sep = "")

fit_dr_nonp_wclsmod <- DR_two_step_wcls_eq2(
  dta = dataset,
  id_varname = "ID",
  decision_time_varname = "Day",
  treatment_varname = "treatment",
  outcome_varname = "number_of_screen1",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL, 
  rand_prob_tilde = 0.6,        
  estimator_initial_value = NULL
)
fit_dr_nonp_wclsmod$beta_hat # 1.49933525   -0.02946211
fit_dr_nonp_wclsmod$beta_se # 0.153213587  0.008847506

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)   
ci_dr_nonp_wclsmod <- round(c(rbind(fit_dr_nonp_wclsmod$beta_hat - t_quantile * fit_dr_nonp_wclsmod$beta_se,
                                    fit_dr_nonp_wclsmod$beta_hat + t_quantile * fit_dr_nonp_wclsmod$beta_se)),3)
# (1.197988,  1.800682)   (-0.04686377, -0.01206044)
p_dr_nonp_wclsmod <- 2 * pnorm(abs(fit_dr_nonp_wclsmod$beta_hat) / fit_dr_nonp_wclsmod$beta_se,  lower.tail = FALSE)
# 3.998510e-20   9.620942e-04

ci_dr_nonp_wclsmod_use1 <- paste("(", paste(ci_dr_nonp_wclsmod[1], ci_dr_nonp_wclsmod[2], sep = ","), ")", sep = "")
ci_dr_nonp_wclsmod_use2 <- paste("(", paste(ci_dr_nonp_wclsmod[3], ci_dr_nonp_wclsmod[4], sep = ","), ")", sep = "")



output1 <- data.frame("Estimator" = c("EMEE", "EMEE-NonP"),
                      "Estimate" = c(fit_wcls$beta_hat, fit_nonp_wcls$beta_hat),
                      "SE" = c(fit_wcls$beta_se_adjusted, fit_nonp_wcls$beta_se_adjusted),
                      "95% CI" = c(ci_wcls_use, ci_nonp_wcls_use),
                      "p-value" = c(p_wcls, p_nonp_wcls))

output2 <- data.frame("Estimator" = c("EMEE", "EMEE-NonP"),
                      "Estimate1" = c(fit_wcls$beta_hat[1], fit_nonp_wcls$beta_hat[1]),
                      "SE1" = c(fit_wcls$beta_se_adjusted[1], fit_nonp_wcls$beta_se_adjusted[1]),
                      "95% CI1" = c(ci_wcls_use1, ci_nonp_wcls_use1),
                      "p-value1" = c(p_wcls[1], p_nonp_wcls[1]),
                      "Estimate2" = c(fit_wcls$beta_hat[2], fit_nonp_wcls$beta_hat[2]),
                      "SE2" = c(fit_wcls$beta_se_adjusted[2], fit_nonp_wcls$beta_se_adjusted[2]),
                      "95% CI2" = c(ci_wcls_use2, ci_nonp_wcls_use2),
                      "p-value2" = c(p_wcls[2], p_nonp_wcls[2]))

output3 <- data.frame("Estimator" = c("EMEE", "EMEE-NonP"),
                      "Estimate1" = c(fit_wclsmod$beta_hat[1], fit_nonp_wclsmod$beta_hat[1]),
                      "SE1" = c(fit_wclsmod$beta_se_adjusted[1], fit_nonp_wclsmod$beta_se_adjusted[1]),
                      "95% CI1" = c(ci_wclsmod_use1, ci_nonp_wclsmod_use1),
                      "p-value1" = c(p_wclsmod[1], p_nonp_wclsmod[1]),
                      "Estimate2" = c(fit_wclsmod$beta_hat[2], fit_nonp_wclsmod$beta_hat[2]),
                      "SE2" = c(fit_wclsmod$beta_se_adjusted[2], fit_nonp_wclsmod$beta_se_adjusted[2]),
                      "95% CI2" = c(ci_wclsmod_use2, ci_nonp_wclsmod_use2),
                      "p-value2" = c(p_wclsmod[2], p_nonp_wclsmod[2]))

library(kableExtra)

# sink("table_generation/simulation_1.txt", append=FALSE)
mycaption <- "Marginal excursion effects of providing push notifications in the Drink Less micro-randomized trial"
latex_code <- kable(output1, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
  # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
  # column_spec(1, bold=T) %>%
  collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)
# sink()

mycaption <- "Marginal excursion effects of providing standard notifications and new notifications in the Drink Less micro-randomized trial"
latex_code <- kable(output2, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
  # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
  # column_spec(1, bold=T) %>%
  collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)

mycaption <- "Moderation effects of the Drink Less micro-randomized trial"
latex_code <- kable(output3, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
  # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
  # column_spec(1, bold=T) %>%
  collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)
