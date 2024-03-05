#data analysis of drink less MRT (count outcome)
library(MASS)
library(parallel)
source("estimator.R")
source("estimator_three_cate.R")
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

EMEE <- fit_EMEE(
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


# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_EMEE <- round(c(rbind(EMEE$beta_hat - t_quantile * EMEE$beta_se,
      EMEE$beta_hat + t_quantile * EMEE$beta_se)),3)
# (0.8593891,  1.3546449)
p_EMEE <- 2 * pt(abs(EMEE$beta_hat) / EMEE$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_EMEE_use <- paste("(", paste(ci_EMEE[1],ci_EMEE[2], sep = ","), ")", sep = "")


EMEE_NonP <- fir_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_EMEE_NonP <- round(c(rbind(EMEE_NonP$beta_hat - t_quantile * EMEE_NonP$beta_se,
      EMEE_NonP$beta_hat + t_quantile * EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_EMEE_NonP <- 2 * pt(abs(EMEE_NonP$beta_hat) / EMEE_NonP$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_EMEE_NonP_use <- paste("(", paste(ci_EMEE_NonP[1], ci_EMEE_NonP[2], sep = ","), ")", sep = "")



DR_EMEE_NonP <- fit_DR_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_DR_EMEE_NonP <- round(c(rbind(DR_EMEE_NonP$beta_hat - t_quantile * DR_EMEE_NonP$beta_se,
                                 DR_EMEE_NonP$beta_hat + t_quantile * DR_EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_DR_EMEE_NonP <- 2 * pt(abs(DR_EMEE_NonP$beta_hat) / DR_EMEE_NonP$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_DR_EMEE_NonP_use <- paste("(", paste(ci_DR_EMEE_NonP[1],ci_DR_EMEE_NonP[2], sep = ","), ")", sep = "")



##### Primary: three treatments, S_t=1#####
prob_A1 <- prob_A2 <- rep(0.3, nrow(dataset))
dataset$prob_A1 <- prob_A1
dataset$prob_A2 <- prob_A2

EMEE <- fit_EMEE(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_EMEE <- round(c(rbind(EMEE$beta_hat - t_quantile * EMEE$beta_se,
                         EMEE$beta_hat + t_quantile * EMEE$beta_se)),3)
ci# (0.8593891,  1.3546449)
p_EMEE <- 2 * pt(abs(fit_EMEE$beta_hat) / fit_EMEE$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_EMEE_use1 <- paste("(", paste(ci_EMEE[1], ci_EMEE[2], sep = ","), ")", sep = "")
ci_EMEE_use2 <- paste("(", paste(ci_EMEE[3], ci_EMEE[4], sep = ","), ")", sep = "")



EMEE_NonP <- fit_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_EMEE_NonP <- round(c(rbind(EMEE_NonP$beta_hat - t_quantile * EMEE_NonP$beta_se,
                              EMEE_NonP$beta_hat + t_quantile * EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_EMEE_NonP <- 2 * pt(abs(EMEE_NonP$beta_hat) / EMEE_NonP$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27


ci_EMEE_NonP_use1 <- paste("(", paste(ci_EMEE_NonP[1], ci_EMEE_NonP[2], sep = ","), ")", sep = "")
ci_EMEE_NonP_use2 <- paste("(", paste(ci_EMEE_NonP[3], ci_EMEE_NonP[4], sep = ","), ")", sep = "")



DR_EMEE_NonP <- fit_DR_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_DR_EMEE_NonP <- round(c(rbind(DR_EMEE_NonP$beta_hat - t_quantile * DR_EMEE_NonP$beta_se,
                              DR_EMEE_NonP$beta_hat + t_quantile * DR_EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_DR_EMEE_NonP <- 2 * pnorm((DR_EMEE_NonP$beta_hat) / DR_EMEE_NonP$beta_se,  lower.tail = FALSE)

ci_DR_EMEE_NonP_use1 <- paste("(", paste(ci_DR_EMEE_NonP[1], ci_DR_EMEE_NonP[2], sep = ","), ")", sep = "")
ci_DR_EMEE_NonP_use2 <- paste("(", paste(ci_DR_EMEE_NonP[3], ci_DR_EMEE_NonP[4], sep = ","), ")", sep = "")


##### Secondary:  S_t = 1 + Day#####
source("estimator.R")
control_vars = c('age','AUDIT_score','Day','gender','Manual_employ', 'Other_employ', 'yesterday_screen_view')
moderator_vars = 'Day'

EMEE <- fit_EMEE(
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


# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_EMEE <- round(c(rbind(EMEE$beta_hat - t_quantile * EMEE$beta_se,
      EMEE$beta_hat + t_quantile * EMEE$beta_se)),3)
# (0.8593891,  1.3546449)
p_EMEE <- 2 * pt(abs(EMEE$beta_hat) / EMEE$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_EMEE_use <- paste("(", paste(ci_EMEE[1],ci_EMEE[2], sep = ","), ")", sep = "")


EMEE_NonP <- fir_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_EMEE_NonP <- round(c(rbind(EMEE_NonP$beta_hat - t_quantile * EMEE_NonP$beta_se,
      EMEE_NonP$beta_hat + t_quantile * EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_EMEE_NonP <- 2 * pt(abs(EMEE_NonP$beta_hat) / EMEE_NonP$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_EMEE_NonP_use <- paste("(", paste(ci_EMEE_NonP[1], ci_EMEE_NonP[2], sep = ","), ")", sep = "")



DR_EMEE_NonP <- fit_DR_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_DR_EMEE_NonP <- round(c(rbind(DR_EMEE_NonP$beta_hat - t_quantile * DR_EMEE_NonP$beta_se,
                                 DR_EMEE_NonP$beta_hat + t_quantile * DR_EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_DR_EMEE_NonP <- 2 * pt(abs(DR_EMEE_NonP$beta_hat) / DR_EMEE_NonP$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_DR_EMEE_NonP_use <- paste("(", paste(ci_DR_EMEE_NonP[1],ci_DR_EMEE_NonP[2], sep = ","), ")", sep = "")



#####S_t = 1 + Yesterday#####
control_vars = c('age','AUDIT_score','Day','gender','Manual_employ', 'Other_employ', 'yesterday_screen_view')
moderator_vars = 'yesterday_screen_view'

EMEE <- fit_EMEE(
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


# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_EMEE <- round(c(rbind(EMEE$beta_hat - t_quantile * EMEE$beta_se,
      EMEE$beta_hat + t_quantile * EMEE$beta_se)),3)
# (0.8593891,  1.3546449)
p_EMEE <- 2 * pt(abs(EMEE$beta_hat) / EMEE$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 6.977297e-17 

ci_EMEE_use <- paste("(", paste(ci_EMEE[1],ci_EMEE[2], sep = ","), ")", sep = "")


EMEE_NonP <- fir_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2)
ci_EMEE_NonP <- round(c(rbind(EMEE_NonP$beta_hat - t_quantile * EMEE_NonP$beta_se,
      EMEE_NonP$beta_hat + t_quantile * EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_EMEE_NonP <- 2 * pt(abs(EMEE_NonP$beta_hat) / EMEE_NonP$beta_se, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_EMEE_NonP_use <- paste("(", paste(ci_EMEE_NonP[1], ci_EMEE_NonP[2], sep = ","), ")", sep = "")



DR_EMEE_NonP <- fit_DR_EMEE_NonP(
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

# estimator, SE, 95% CI, p-value
t_quantile <- qnorm(1 - 0.05/2) 
ci_DR_EMEE_NonP <- round(c(rbind(DR_EMEE_NonP$beta_hat - t_quantile * DR_EMEE_NonP$beta_se,
                                 DR_EMEE_NonP$beta_hat + t_quantile * DR_EMEE_NonP$beta_se)), 3)
# (0.9273521,  1.3002999)
p_DR_EMEE_NonP <- 2 * pt(abs(DR_EMEE_NonP$beta_hat) / DR_EMEE_NonP$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 4.791444e-27

ci_DR_EMEE_NonP_use <- paste("(", paste(ci_DR_EMEE_NonP[1],ci_DR_EMEE_NonP[2], sep = ","), ")", sep = "")
