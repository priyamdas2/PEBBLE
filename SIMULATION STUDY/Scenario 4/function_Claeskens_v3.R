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
# ########### Generating Data ##################################
# set.seed(1)
# Num_of_exp <- 10
# n <- 100
# p <- 5
# Beta <- 10*(runif(p)-0.5)
# 
# 
# ##### BOOTSTRAP parameters ##################################
# c <- 1.5
# beta_a <- 0.5
# beta_b <- 1.5
# coord_num <- 2
# num_BOOT <- 1000
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
# 
# #############################################################
# 
# n <- dim(X)[1]
# p <- dim(X)[2]
# temp_prod <- exp(X%*%Beta_hat)
# Psi <- diag(as.vector((t(t(Y)) - temp_prod/((1+temp_prod)))))%*%X
# temp_factor_array <- -temp_prod/((1+temp_prod)^2)
# Psi_dot <- array(NA,dim=c(p,p,n))
# Xi <- temp_factor_array
# Xi_dot <- matrix(0,n,p)
# 
# for(i in 1:n)
# {Psi_dot[,,i] <- temp_factor_array[i]*X[i,]%*%t(X[i,])
#  temp_quan <- -temp_prod[i]*(1-temp_prod[i])/((1+temp_prod[i])^3)
#  Xi_dot[i,] <- temp_quan*X[i,]
# }
# SRSWR_index <- rdunif(n,1,n)
# 




T_star_Claeskens <- function(X,Y,X_0,Beta_hat,SRSWR_index,Psi,Psi_dot,Xi_dot)
{n <- dim(X)[1]
p <- dim(X)[2]

sum_Psi_dot_star <- matrix(0,p,p)
sum_Psi_star <- matrix(0,p,1)


for(i in 1:n)
{sum_Psi_star <- sum_Psi_star + t(t(Psi[SRSWR_index[i],]))
sum_Psi_dot_star <- sum_Psi_dot_star + Psi_dot[,,SRSWR_index[i]]
}
sum_Psi_dot_star_temp_neg <- nearPD(-sum_Psi_dot_star)
inv_sum_Psi_dot_star_neg <- solve(sum_Psi_dot_star_temp_neg$mat)
U_star <- inv_sum_Psi_dot_star_neg%*%sum_Psi_star
Beta_3star <- t(t(Beta_hat)) + U_star
T_3star <- sqrt(n)*(exp(sum(X_0*Beta_3star)) - exp(sum(X_0*Beta_hat)))


######### H_4star Calculation #################
Xi_dot_star <- Xi_dot[SRSWR_index,]
X_star <- X[SRSWR_index,]
right_factor <- matrix(0,n,p)


right_factor <- matrix(0,p,1)
for(i in 1:n)
{right_factor <- right_factor + sum(X_star[i,]*U_star)*t(t(Xi_dot_star[i,]))}
Beta_4star <- Beta_3star + 0.5*inv_sum_Psi_dot_star_neg%*%right_factor
T_4star <- sqrt(n)*(exp(sum(X_0*Beta_4star)) - exp(sum(X_0*Beta_hat)))

return(list(T_3star = T_3star, T_4star = T_4star))
}


H_star_Claeskens <- function(X,Y,Beta_hat,SRSWR_index,Psi,Psi_dot,Xi_dot)
{n <- dim(X)[1]
p <- dim(X)[2]

sum_Psi_dot_star <- matrix(0,p,p)
sum_Psi_star <- matrix(0,p,1)


for(i in 1:n)
{sum_Psi_star <- sum_Psi_star + t(t(Psi[SRSWR_index[i],]))
sum_Psi_dot_star <- sum_Psi_dot_star + Psi_dot[,,SRSWR_index[i]]
}
sum_Psi_dot_star_temp_neg <- nearPD(-sum_Psi_dot_star)
inv_sum_Psi_dot_star_neg <- solve(sum_Psi_dot_star_temp_neg$mat)
U_star <- inv_sum_Psi_dot_star_neg%*%sum_Psi_star
Beta_3star <- t(t(Beta_hat)) + U_star
H_3star <- sqrt(n)*(Beta_3star - t(t(Beta_hat)))


######### H_4star Calculation #################
Xi_dot_star <- Xi_dot[SRSWR_index,]
X_star <- X[SRSWR_index,]
right_factor <- matrix(0,n,p)


right_factor <- matrix(0,p,1)
for(i in 1:n)
{right_factor <- right_factor + sum(X_star[i,]*U_star)*t(t(Xi_dot_star[i,]))}
Beta_4star <- Beta_3star + 0.5*inv_sum_Psi_dot_star_neg%*%right_factor
H_4star <- sqrt(n)*(Beta_4star - t(t(Beta_hat)))

return(list(H_3star = H_3star, H_4star = H_4star))
}

#############################################
#############################################
#############################################

bounds_Claeskens <- function(X,Y,X_0,num_BOOT,x0,iter_no)
{n <- dim(X)[1]
p <- dim(X)[2]

#### Finding Beta_hat and L_hat ##############
negloglikelogistic_ORIGINAL_BETA <- function(Beta){
  return(-loglikelogistic_ORIGINAL(Beta, Y, X))
}
x0 <- array(1,p)

dd <- optimr(x0,negloglikelogistic_ORIGINAL_BETA)
Beta_hat <- dd$par

#######Finding Psi,Psi_dot,Xi,Xi_dot #########

n <- dim(X)[1]
p <- dim(X)[2]
temp_prod <- exp(X%*%Beta_hat)
Psi <- diag(as.vector((t(t(Y)) - temp_prod/((1+temp_prod)))))%*%X
temp_factor_array <- -temp_prod/((1+temp_prod)^2)
Psi_dot <- array(NA,dim=c(p,p,n))
Xi <- temp_factor_array
Xi_dot <- matrix(0,n,p)

for(i in 1:n)
{Psi_dot[,,i] <- temp_factor_array[i]*X[i,]%*%t(X[i,])
temp_quan <- -temp_prod[i]*(1-temp_prod[i])/((1+temp_prod[i])^3)
Xi_dot[i,] <- temp_quan*X[i,]
}
SRSWR_index <- rdunif(n,1,n)
#############################################
H_3star_MATRIX <- matrix(NA,num_BOOT,p)
H_4star_MATRIX <- matrix(NA,num_BOOT,p)

H_3star_norms <- array(NA,num_BOOT)
H_4star_norms <- array(NA,num_BOOT)

T_3star_array <- array(NA,num_BOOT)
T_4star_array <- array(NA,num_BOOT)

for(i in 1:num_BOOT)
{ #print(i)
  # message(paste('Claeskens Bootstrap iteration number =',i,',iter_no = ',iter_no))
  SRSWR_index <- rdunif(n,1,n)
  DD <- H_star_Claeskens(X,Y,Beta_hat,SRSWR_index,Psi,Psi_dot,Xi_dot)
  H_3star_MATRIX[i,] <- as.array(DD$H_3star)
  H_4star_MATRIX[i,] <- as.array(DD$H_4star)
  H_3star_norms[i] <- sum(as.array(DD$H_3star)^2)
  H_4star_norms[i] <- sum(as.array(DD$H_4star)^2)
  
  DDDD <- T_star_Claeskens(X,Y,X_0,Beta_hat,SRSWR_index,Psi,Psi_dot,Xi_dot)
  T_3star_array[i] <- DDDD$T_3star
  T_4star_array[i] <- DDDD$T_4star 
}

index_q_05 <- array(NA,p)
index_q_10 <- array(NA,p)
index_q_90 <- array(NA,p)
index_q_95 <- array(NA,p)

H_3star_j_05 <- array(NA,p)
H_3star_j_10 <- array(NA,p)
H_3star_j_90 <- array(NA,p)
H_3star_j_95 <- array(NA,p)


index_q_05_alt <- array(NA,p)
index_q_10_alt <- array(NA,p)
index_q_90_alt <- array(NA,p)
index_q_95_alt <- array(NA,p)

H_4star_j_05 <- array(NA,p)
H_4star_j_10 <- array(NA,p)
H_4star_j_90 <- array(NA,p)
H_4star_j_95 <- array(NA,p)

temp_common <- array(NA,p)

for(j in 1:p)
{index_q_05[j] <- which.min(abs(H_3star_MATRIX[,j] - quantile(H_3star_MATRIX[,j], .05)))
index_q_10[j] <- which.min(abs(H_3star_MATRIX[,j] - quantile(H_3star_MATRIX[,j], .1)))
index_q_90[j] <- which.min(abs(H_3star_MATRIX[,j] - quantile(H_3star_MATRIX[,j], .90)))
index_q_95[j] <- which.min(abs(H_3star_MATRIX[,j] - quantile(H_3star_MATRIX[,j], .95)))

H_3star_j_05[j] <- H_3star_MATRIX[index_q_05[j],j]
H_3star_j_10[j] <- H_3star_MATRIX[index_q_10[j],j]
H_3star_j_90[j] <- H_3star_MATRIX[index_q_90[j],j]
H_3star_j_95[j] <- H_3star_MATRIX[index_q_95[j],j]

index_q_05_alt[j] <- which.min(abs(H_4star_MATRIX[,j] - quantile(H_4star_MATRIX[,j], .05)))
index_q_10_alt[j] <- which.min(abs(H_4star_MATRIX[,j] - quantile(H_4star_MATRIX[,j], .1)))
index_q_90_alt[j] <- which.min(abs(H_4star_MATRIX[,j] - quantile(H_4star_MATRIX[,j], .90)))
index_q_95_alt[j] <- which.min(abs(H_4star_MATRIX[,j] - quantile(H_4star_MATRIX[,j], .95)))

H_4star_j_05[j] <- H_4star_MATRIX[index_q_05_alt[j],j]
H_4star_j_10[j] <- H_4star_MATRIX[index_q_10_alt[j],j]
H_4star_j_90[j] <- H_4star_MATRIX[index_q_90_alt[j],j]
H_4star_j_95[j] <- H_4star_MATRIX[index_q_95_alt[j],j]

}


T_3star_j_05 <- quantile(T_3star_array,.05)
T_3star_j_10 <- quantile(T_3star_array,.10)
T_3star_j_90 <- quantile(T_3star_array,.90)
T_3star_j_95 <- quantile(T_3star_array,.95)

T_3star_quantiles <- c(T_3star_j_05,T_3star_j_10,T_3star_j_90,T_3star_j_95)

T_4star_j_05 <- quantile(T_4star_array,.05)
T_4star_j_10 <- quantile(T_4star_array,.10)
T_4star_j_90 <- quantile(T_4star_array,.90)
T_4star_j_95 <- quantile(T_4star_array,.95)

T_4star_quantiles <- c(T_4star_j_05,T_4star_j_10,T_4star_j_90,T_4star_j_95)


upper_bound <- Beta_hat - H_3star_j_05/(sqrt(n))
lower_bound <- Beta_hat - H_3star_j_95/(sqrt(n))
width <- upper_bound - lower_bound
upper_bound_90 <- Beta_hat - H_3star_j_10/(sqrt(n))
lower_bound_10 <- Beta_hat - H_3star_j_90/(sqrt(n))


upper_bound_alt <- Beta_hat - H_4star_j_05/(sqrt(n))
lower_bound_alt <- Beta_hat - H_4star_j_95/(sqrt(n))
width_alt <- upper_bound_alt - lower_bound_alt
upper_bound_90_alt <- Beta_hat - H_4star_j_10/(sqrt(n))
lower_bound_10_alt <- Beta_hat - H_4star_j_90/(sqrt(n))

H_3star_whole_05_10_90_95 <- array(NA,4)
H_4star_whole_05_10_90_95 <- array(NA,4)

index_whole <- which.min(abs(H_3star_norms - quantile(H_3star_norms, .05)))
index_whole[2] <- which.min(abs(H_3star_norms - quantile(H_3star_norms, .1)))
index_whole[3] <- which.min(abs(H_3star_norms - quantile(H_3star_norms, .9)))
index_whole[4] <- which.min(abs(H_3star_norms - quantile(H_3star_norms, .95)))

H_3star_whole_05_10_90_95 <- H_3star_norms[index_whole]

index_whole_alt <- which.min(abs(H_4star_norms - quantile(H_4star_norms, .05)))
index_whole_alt[2] <- which.min(abs(H_4star_norms - quantile(H_4star_norms, .1)))
index_whole_alt[3] <- which.min(abs(H_4star_norms - quantile(H_4star_norms, .9)))
index_whole_alt[4] <- which.min(abs(H_4star_norms - quantile(H_4star_norms, .95)))

H_4star_whole_05_10_90_95 <- H_4star_norms[index_whole_alt]



return(list(lower_bound = lower_bound, upper_bound = upper_bound,lower_bound_10 = lower_bound_10,
            upper_bound_90 = upper_bound_90, width = width,H_3star_whole_05_10_90_95 = H_3star_whole_05_10_90_95,
            lower_bound_alt = lower_bound_alt, upper_bound_alt = upper_bound_alt,lower_bound_10_alt = lower_bound_10_alt,
            upper_bound_90_alt = upper_bound_90_alt, width_alt = width_alt, H_4star_whole_05_10_90_95 = H_4star_whole_05_10_90_95,
            T_3star_quantiles = T_3star_quantiles, T_4star_quantiles = T_4star_quantiles))

}

# bounds_Claeskens(X,Y,num_BOOT,x0)