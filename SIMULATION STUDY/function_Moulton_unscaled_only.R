# rm(list=ls())
# library(Rcpp)
# library(RcppArmadillo)
# library(DEoptim)
# library(optimr)
# library(expm)
# library(Rfast)
# library(purrr)
# 
# setwd("D:/Logistic regression Project")
# sourceCpp("function_supp.cpp")
# sourceCpp("function_rough.cpp")
# source("function_supp_R.R")
# 
# 
# 
# ############ Generating Data ##################################
# set.seed(2)
# Num_of_exp <- 2
# n <- 200
# p <- 3
# Beta <- 3*(runif(p)-0.5)
# 
# 
# ##### BOOTSTRAP parameters ##################################
# c <- 1.5
# beta_a <- 0.5
# beta_b <- 1.5
# coord_num <- 2
# num_BOOT <- 3
# #############################################################
# 
# 
# X <- matrix(runif(n*p)-0.5,n,p)
# P_x_true <- exp(X%*%Beta)/(1+exp(X%*%Beta))
# Y <- array(NA,n)
# for(ii in 1:n)
# {Y[ii] <- rbinom(1,1,P_x_true[ii])}
# 
# 
# ### Regular logistic regression ###################
# negloglikelogistic_ORIGINAL_BETA <- function(Beta){
#   return(-loglikelogistic_ORIGINAL(Beta, Y, X))
# }
# 
# x0 <- array(1,p)
# dd <- optimr(x0,negloglikelogistic_ORIGINAL_BETA)
# Beta_hat <- dd$par

#############################################################


H_star_Moulton <- function(X,Y,Beta_hat,SRSWR_index, e_array, G_T_G_inv_G_T)
{n <- dim(X)[1]
p <- dim(X)[2]
e_star_array <- e_array[SRSWR_index]
Beta_star_star <- t(t(Beta_hat))+(1/n)*G_T_G_inv_G_T%*%t(t(e_star_array))

H_star_Moulton_unscaled <- sqrt(n)*(Beta_star_star - t(t(Beta_hat)))

return(H_star_Moulton_unscaled)
}

########################################################
########################################################
########################################################

bounds_Moulton <- function(X,Y,num_BOOT,x0,iter_no)
{n <- dim(X)[1]
p <- dim(X)[2]

#### Finding Beta_hat and L_hat ##############
negloglikelogistic_ORIGINAL_BETA <- function(Beta){
  return(-loglikelogistic_ORIGINAL(Beta, Y, X))
}
x0 <- array(1,p)

dd <- optimr(x0,negloglikelogistic_ORIGINAL_BETA)
Beta_hat <- dd$par

###############################################

temp_prod <- exp(X%*%Beta_hat)
s <- t(t(Y)) - temp_prod/((1+temp_prod))
V <- temp_prod/((1+temp_prod)^2)
V_pow_half <- sqrtm(diag(as.vector(V)))
G <- V_pow_half%*%X

### Finding L_hat ##########
coef_i <- temp_prod/((1+temp_prod)^2)
temp_mat_sum <- matrix(0, p, p)
for(i in 1:n)
{temp_mat_sum <- temp_mat_sum + coef_i[i]*t(t(X[i,]))%*%X[i,]}

L_hat <- temp_mat_sum/n
L_hat_obj <- nearPD(L_hat)
L_hat <- L_hat_obj$mat
L_hat_inv <- solve(L_hat)
#############################

h_mat <- (G%*%(L_hat_inv)%*%t(G))/n
e_array <- array(NA,n)
for(i in 1:n)
{e_array[i] <- s[i]/(abs(V[i]*(1-h_mat[i,i]))^0.5)}
G_T_G_inv_G_T <- L_hat_inv%*%t(G)


###############################################


H_star_Moulton_alt_MATRIX <- matrix(NA,num_BOOT,p)
H_star_Moulton_alt_norms <- array(NA,num_BOOT)

for(i in 1:num_BOOT)
{ #print(i)
  #message(paste('Moulton Bootstrap iteration number =',i,',iter_no = ',iter_no))
  SRSWR_index <- rdunif(n,1,n)
  DD <- H_star_Moulton(X,Y,Beta_hat,SRSWR_index, e_array, G_T_G_inv_G_T)
  H_star_Moulton_alt_MATRIX[i,] <- as.array(DD)
  
  H_star_Moulton_alt_norms[i] <- sum(as.array(DD)^2)
}


index_q_05_alt <- array(NA,p)
index_q_10_alt <- array(NA,p)
index_q_90_alt <- array(NA,p)
index_q_95_alt <- array(NA,p)

H_star_j_05_alt <- array(NA,p)
H_star_j_10_alt <- array(NA,p)
H_star_j_90_alt <- array(NA,p)
H_star_j_95_alt <- array(NA,p)

temp_common <- array(NA,p)

for(j in 1:p)
{index_q_05_alt[j] <- which.min(abs(H_star_Moulton_alt_MATRIX[,j] - quantile(H_star_Moulton_alt_MATRIX[,j], .05)))
index_q_10_alt[j] <- which.min(abs(H_star_Moulton_alt_MATRIX[,j] - quantile(H_star_Moulton_alt_MATRIX[,j], .1)))
index_q_90_alt[j] <- which.min(abs(H_star_Moulton_alt_MATRIX[,j] - quantile(H_star_Moulton_alt_MATRIX[,j], .90)))
index_q_95_alt[j] <- which.min(abs(H_star_Moulton_alt_MATRIX[,j] - quantile(H_star_Moulton_alt_MATRIX[,j], .95)))

H_star_j_05_alt[j] <- H_star_Moulton_alt_MATRIX[index_q_05_alt[j],j]
H_star_j_10_alt[j] <- H_star_Moulton_alt_MATRIX[index_q_10_alt[j],j]
H_star_j_90_alt[j] <- H_star_Moulton_alt_MATRIX[index_q_90_alt[j],j]
H_star_j_95_alt[j] <- H_star_Moulton_alt_MATRIX[index_q_95_alt[j],j]
}


upper_bound_alt <- Beta_hat - H_star_j_05_alt/(sqrt(n))
lower_bound_alt <- Beta_hat - H_star_j_95_alt/(sqrt(n))
width_alt <- upper_bound_alt - lower_bound_alt
upper_bound_90_alt <- Beta_hat - H_star_j_10_alt/(sqrt(n))
lower_bound_10_alt <- Beta_hat - H_star_j_90_alt/(sqrt(n))



H_star_Moulton_alt_whole_05_10_90_95 <- array(NA,4)

index_whole <- which.min(abs(H_star_Moulton_alt_norms - quantile(H_star_Moulton_alt_norms, .05)))
index_whole[2] <- which.min(abs(H_star_Moulton_alt_norms - quantile(H_star_Moulton_alt_norms, .1)))
index_whole[3] <- which.min(abs(H_star_Moulton_alt_norms - quantile(H_star_Moulton_alt_norms, .9)))
index_whole[4] <- which.min(abs(H_star_Moulton_alt_norms - quantile(H_star_Moulton_alt_norms, .95)))

H_star_Moulton_alt_whole_05_10_90_95 <- H_star_Moulton_alt_norms[index_whole]


return(list(lower_bound_alt = lower_bound_alt, upper_bound_alt = upper_bound_alt,lower_bound_10_alt = lower_bound_10_alt,
            upper_bound_90_alt = upper_bound_90_alt, width_alt = width_alt,
            H_star_Moulton_alt_whole_05_10_90_95 = H_star_Moulton_alt_whole_05_10_90_95))

}
