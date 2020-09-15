rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(DEoptim)
library(optimr)
library(expm)
library(Rfast)
library(purrr)
library(MASS)


setwd("D:/Logistic regression Project/GITHUB CODES/REAL DATA ANALYSIS")
sourceCpp("function_supp.cpp")
sourceCpp("function_rough.cpp")
source("function_Pebble_REAL_DATA.R")


ptm <- proc.time()

############ Generating Data ##################################
seed_no <- 2
set.seed(seed_no)


MAT <- as.matrix(read.table("Real_data_reformat.csv", header = TRUE, sep=","))

X_temp <- MAT[,1:5]


Y <- MAT[,6]

Num_of_exp <- 1
X <- cbind(rep(1,dim(X_temp)[2]), X_temp)

n <- dim(X)[1]
p <- dim(X)[2]

##### BOOTSTRAP parameters ##################################
c <- 1.0
beta_a <- 0.5
beta_b <- 1.5
num_BOOT <- 1000

p_1 <- max(p + 1,4)
b <- sqrt(c*n^(-1/p_1))
# coord_num <- 2
bound_max <- 10000

SIGMA_BOUND <- 50
Norm_var <- 0.25
#############################################################




Beta_hat_MAT <- matrix(NA, Num_of_exp,p)
Beta_star_MAT <- matrix(NA, Num_of_exp,p)

Width_Pebble <- matrix(NA,Num_of_exp,p)

H <- matrix(NA,Num_of_exp,p)
H_norm <- array(NA,Num_of_exp)

z_norm_star_record_ALL <- list()
z_norm_ALL <- matrix(NA,Num_of_exp,p)



iter_no = 1
print(iter_no)
  set.seed(iter_no)
  
  ### Regular logistic regression ###################
  negloglikelogistic_ORIGINAL_BETA <- function(Beta){
    return(-loglikelogistic_ORIGINAL(Beta, Y, X))
  }
  
  x0 <- array(0,p)
  dd <- optimr(x0,negloglikelogistic_ORIGINAL_BETA)
  Beta_hat <- dd$par
  Beta_hat_MAT[iter_no,] <- t(Beta_hat)
  
  ################### FINDING H (true) ###################
  
  temp_prod <- exp(X%*%Beta_hat)
  coef_i <- temp_prod/((1+temp_prod)^2)  # For hat(L_n)
  
  temp_prod_2 <- temp_prod
  p_hat <- temp_prod_2/(1+temp_prod_2)  # For hat(M_n)
  
  
  temp_mat_sum <- matrix(0, p, p)
  temp_mat_sum_2 <- matrix(0, p, p)
  
  for(i in 1:n)
  {prod_temp <- t(t(X[i,]))%*%X[i,]
  
  temp_mat_sum <- temp_mat_sum + prod_temp*coef_i[i]   # For hat(L_n)
  temp_mat_sum_2 <- temp_mat_sum_2 + prod_temp*(Y[i] - p_hat[i])^2}   # For hat(M_n)
  
  # \hat{\Sigma} = \hat{L}^{-1}\hat{M}\hat{L}^{-1}
  
  L_hat <- temp_mat_sum/n
  #L_hat <- L_hat + diag(p)/sqrt(n)
  L_hat_obj <- nearPD(L_hat)
  L_hat <- L_hat_obj$mat
  L_hat_inv <- solve(L_hat)
  
  M_hat <- temp_mat_sum_2/n
  M_hat_obj <- nearPD(M_hat)
  M_hat <- M_hat_obj$mat
  
  Sigma_hat <- L_hat_inv%*%M_hat%*%L_hat_inv
  
  # \hat{\Sigma}_{j,n} = (j,j) th elemet of \hat{\Sigma}
  
  Sigma_hat_obj <- nearPD(Sigma_hat)
  Sigma_hat <- Sigma_hat_obj$mat
  Sigma_hat_inv <- solve(Sigma_hat)
  Sigma_hat_pow_neg_half <- sqrtm(Sigma_hat_inv)
  Sigma_hat_diag <- diag(Sigma_hat)
  #########################################################
  
  z_norm <- t(rnorm(p,0,Norm_var))
  
  #########################################################
  ###### Pebble Bootstrap logistic regression ############
  #########################################################
  alpha <- 0.1
  DDD <- bounds_Pebble_alpha(X,Y,b, num_BOOT, beta_a, beta_b,x0,iter_no,Sigma_hat_diag,Norm_var,alpha)
  H_star_j_alpha_by_2 <- DDD$H_star_j_alpha_by_2
  H_star_j_alpha <- DDD$H_star_j_alpha
  H_star_j_one_minus_alpha <- DDD$H_star_j_one_minus_alpha
  H_star_j_one_minus_alpha_by_2 <- DDD$H_star_j_one_minus_alpha_by_2

  
  l_star <- rep(NA,p)
  u_star <- rep(NA,p)
  lb <- rep(NA,p)
  ub <- rep(NA,p)
  
  u_star_for_upper_CI <- rep(NA,p)
  lb_for_upper_CI <- rep(NA,p)
  
  l_star_for_lower_CI <- rep(NA,p)
  ub_for_lower_CI <- rep(NA,p)
  
for(j in 1:p)
{quantity_here <- b*(Sigma_hat_diag[j]^(-1/2))*L_hat_inv[j,]%*%t(z_norm)
 l_star[j] <- H_star_j_alpha_by_2[j] - quantity_here
 u_star[j] <- H_star_j_one_minus_alpha_by_2[j] - quantity_here
 u_star_for_upper_CI[j] <- H_star_j_one_minus_alpha[j] - quantity_here
 l_star_for_lower_CI[j] <- H_star_j_alpha[j] - quantity_here
 lb[j] <- Beta_hat[j] - (Sigma_hat_diag[j]^(1/2))*u_star[j]/sqrt(n)
 ub[j] <- Beta_hat[j] - (Sigma_hat_diag[j]^(1/2))*l_star[j]/sqrt(n)
 lb_for_upper_CI[j] <- Beta_hat[j] - (Sigma_hat_diag[j]^(1/2))*u_star_for_upper_CI[j]/sqrt(n)
 ub_for_lower_CI[j] <- Beta_hat[j] - (Sigma_hat_diag[j]^(1/2))*l_star_for_lower_CI[j]/sqrt(n)
}
  
Beta_hat
lb
ub
lb_for_upper_CI
ub_for_lower_CI
  
summary_mat <- cbind(Beta_hat,lb,ub,lb_for_upper_CI,ub_for_lower_CI)
#########################################################
  





write.csv(round(summary_mat,3), paste0("Real_data_summary.csv"), row.names =  c("intercept",colnames(MAT[,1:5])))

