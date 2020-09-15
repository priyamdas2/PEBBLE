rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(DEoptim)
library(optimr)
library(expm)
library(Rfast)
library(purrr)
library(MASS)


setwd("D:/Logistic regression Project/GITHUB CODES")
sourceCpp("function_supp.cpp")
sourceCpp("function_rough.cpp")
source("function_Pebble.R")
source("function_Moulton_unscaled_only.R")
source("function_Claeskens.R")

ptm <- proc.time()

############ Generating Data ##################################
seed_no <- 2
set.seed(seed_no)
Num_of_exp <- 3
n <- 50                     # Change (n,p) as required
p <- 3


Beta_whole <- c(1,.5,-2,-0.75, 1.5,-1, 1.85, -1.6, -1.75, -0.75)
Beta <- Beta_whole[1:p]

X_SIGMA <- diag(p)
for(i in 1:p)
{for(j in 1:p)
{X_SIGMA[i,j] <- 0.5^abs(i-j)}}

X <- mvrnorm(n, rep(0, p), X_SIGMA)
P_x_true <- exp(X%*%Beta)/(1+exp(X%*%Beta))

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

Coverage_Pebble_middle <- matrix(NA,Num_of_exp,p)
Coverage_Pebble_upper <- matrix(NA,Num_of_exp,p) # Upper 90% CI
Coverage_Pebble_lower <- matrix(NA,Num_of_exp,p) # Lower 90% CI
Coverage_Pebble_Hnorm_mid_up_low <- matrix(NA,Num_of_exp,3)

Coverage_Normal_middle <- matrix(0,Num_of_exp,p)
Coverage_Normal_upper <- matrix(0,Num_of_exp,p) # Upper 90% CI
Coverage_Normal_lower <- matrix(0,Num_of_exp,p) # Lower 90% CI
Coverage_Normal_Hnorm_mid_up_low <- matrix(0,Num_of_exp,3)

Coverage_Moulton_alt_middle <- matrix(0,Num_of_exp,p)
Coverage_Moulton_alt_upper <- matrix(0,Num_of_exp,p) # Upper 90% CI
Coverage_Moulton_alt_lower <- matrix(0,Num_of_exp,p) # Lower 90% CI
Coverage_Moulton_alt_Hnorm_mid_up_low <- matrix(0,Num_of_exp,3)

Coverage_Claeskens_middle <- matrix(0,Num_of_exp,p)
Coverage_Claeskens_upper <- matrix(0,Num_of_exp,p) # Upper 90% CI
Coverage_Claeskens_lower <- matrix(0,Num_of_exp,p) # Lower 90% CI
Coverage_Claeskens_Hnorm_mid_up_low <- matrix(0,Num_of_exp,3)

Coverage_Claeskens_alt_middle <- matrix(0,Num_of_exp,p)
Coverage_Claeskens_alt_upper <- matrix(0,Num_of_exp,p) # Upper 90% CI
Coverage_Claeskens_alt_lower <- matrix(0,Num_of_exp,p) # Lower 90% CI
Coverage_Claeskens_alt_Hnorm_mid_up_low <- matrix(0,Num_of_exp,3)



Width_Pebble <- matrix(NA,Num_of_exp,p)
Width_Normal <- matrix(NA,Num_of_exp,p)
Width_Moulton <- matrix(NA,Num_of_exp,p)
Width_Moulton_alt <- matrix(NA,Num_of_exp,p)
Width_Claeskens <- matrix(NA,Num_of_exp,p)
Width_Claeskens_alt <- matrix(NA,Num_of_exp,p)



H <- matrix(NA,Num_of_exp,p)
H_norm <- array(NA,Num_of_exp)

z_norm_star_record_ALL <- list()
z_norm_ALL <- matrix(NA,Num_of_exp,p)



# iter_no = 1
for(iter_no in 1:Num_of_exp)
{ print("Experiment Number")
  print(iter_no)
  set.seed(iter_no)
  
  fact_check <- 0
  
  while(fact_check == 0)
  {Y <- array(NA,n)
  for(ii in 1:n)
  {Y[ii] <- rbinom(1,1,P_x_true[ii])}
  
  ### Regular logistic regression ###################
  negloglikelogistic_ORIGINAL_BETA <- function(Beta){
    return(-loglikelogistic_ORIGINAL(Beta, Y, X))
  }
  
  x0 <- array(1,p)
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
  
  if(max(Sigma_hat_diag) < SIGMA_BOUND)
  {fact_check <- 1}
  }
  
  
  #########################################################
  
  z_norm <- t(rnorm(p,0,Norm_var))
  
  #########################################################
  ###### Pebble Bootstrap logistic regression ############
  #########################################################
  
  DDD <- bounds_Pebble(X,Y,b, num_BOOT, beta_a, beta_b,x0,iter_no,Sigma_hat_diag,Norm_var)
  H_star_j_95 <- DDD$H_star_j_95
  H_star_j_90 <- DDD$H_star_j_90
  H_star_j_10 <- DDD$H_star_j_10
  H_star_j_05 <- DDD$H_star_j_05
  H_star_whole_05_10_90_95 <- DDD$H_star_whole_05_10_90_95
  z_norm_star_record_ALL[[iter_no]] <- DDD$z_norm_star_record
  #########################################################
  
  if(is.na(H_star_j_95[1]) == FALSE)
  {Temp_array_1 <- array(NA,p)
  
  # Temp_array_1 = (sigma_hat[jjj]^0.5)/sqrt(n)*H_star[j]
  for(jjj in 1:p)
  {Temp_array_1[jjj] <- (Sigma_hat_diag[jjj]^(-1/2))*(sqrt(n)*(Beta_hat[jjj] - Beta[jjj]) + b*(L_hat_inv[jjj,]%*%t(z_norm)))}
  
  while(max(abs(Temp_array_1)) > bound_max)
  {z_norm <- t(rnorm(p))
  for(jjj in 1:p)
  {Temp_array_1[jjj] <- (Sigma_hat_diag[jjj]^(-1/2))*(sqrt(n)*(Beta_hat[jjj] - Beta[jjj]) + b*(L_hat_inv[jjj,]%*%t(z_norm)))}
  }
  
  z_norm_ALL[iter_no,] <- z_norm
  
  #Temp_array_2 <- sqrt(n)*(Sigma_hat_diag^(1/2))*Temp_array_1
  H[iter_no,] <- as.array(Temp_array_1)
  
  H_hat_whole <- Sigma_hat_pow_neg_half%*%(sqrt(n)*t(t(Beta_hat - Beta))+b*L_hat_inv%*%t(z_norm))
  H_hat_norm <- sum(H_hat_whole^2)
  
  #######################################################
  
  ######Coverage of H_n norm ##########################
  
  if(H_hat_norm <= H_star_whole_05_10_90_95[4] && H_hat_norm>= H_star_whole_05_10_90_95[1])
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,1] <- 1}
  else
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,1] <- 0}
  
  if(H_hat_norm >= H_star_whole_05_10_90_95[2])
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,2] <- 1}
  else
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,2] <- 0}
  
  if(H_hat_norm <= H_star_whole_05_10_90_95[3])
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,3] <- 1}
  else
  {Coverage_Pebble_Hnorm_mid_up_low[iter_no,3] <- 0}
  
  ############  DIRECT Coverage using H_n, ###################
  
  ### Middle 90% CI ######
  for(j in 1:p)
  {if(H[iter_no,j] <= H_star_j_95[j] && H[iter_no,j] >= H_star_j_05[j])
  {Coverage_Pebble_middle[iter_no,j] <- 1}
    else{Coverage_Pebble_middle[iter_no,j] <- 0}
    
    
    ##  Upper 90% CI #######
    if(H[iter_no,j] >= H_star_j_10[j])
    {Coverage_Pebble_upper[iter_no,j] <- 1}
    else
    {Coverage_Pebble_upper[iter_no,j] <- 0}
    
    ##  Lower 90% CI #######
    if(H[iter_no,j] <= H_star_j_90[j])
    {Coverage_Pebble_lower[iter_no,j] <- 1}
    else
    {Coverage_Pebble_lower[iter_no,j] <- 0}
  }
  Sigma_hat_pow_half <- sqrtm(Sigma_hat)
  Width_Pebble[iter_no,] <- as.array((Sigma_hat_diag^(1/2))*(H_star_j_95 - H_star_j_05)/sqrt(n))
  }
  
  #########################################################
  ############  Normal Logistic Regression  ###############
  #########################################################
  
  Sigma_hat_here <- L_hat_inv
  Sigma_hat_here_inv <- solve(Sigma_hat_here)
  Sigma_hat_here_pow_neg_half <- sqrtm(Sigma_hat_here_inv)
  diag_Sigma_hat_here <- diag(Sigma_hat_here)
  temp_TRUE <- Beta
  temp_low <- Beta_hat - (diag_Sigma_hat_here^(1/2))*qnorm(0.95)/(sqrt(n))
  temp_up <- Beta_hat - (diag_Sigma_hat_here^(1/2))*qnorm(0.05)/(sqrt(n))
  temp_low_10 <- Beta_hat - (diag_Sigma_hat_here^(1/2))*qnorm(0.9)/(sqrt(n))
  temp_up_90 <- Beta_hat - (diag_Sigma_hat_here^(1/2))*qnorm(0.1)/(sqrt(n))
  
  H_hat_whole_normal <- Sigma_hat_here_pow_neg_half%*%(sqrt(n)*t(t(Beta_hat - Beta)))
  H_hat_norm <- sum(H_hat_whole_normal^2)
  
  q_95 <- qchisq(.95, df=p)
  q_90 <- qchisq(.90, df=p)
  q_10 <- qchisq(.10, df=p)
  q_05 <- qchisq(.05, df=p)
  
  ###### Coverage of H_n norm ##########################
  
  if(H_hat_norm <= q_95 && H_hat_norm>= q_05)
  {Coverage_Normal_Hnorm_mid_up_low[iter_no,1] <- 1}
  
  if(H_hat_norm >= q_10)
  {Coverage_Normal_Hnorm_mid_up_low[iter_no,2] <- 1}
  
  if(H_hat_norm <= q_90)
  {Coverage_Normal_Hnorm_mid_up_low[iter_no,3] <- 1}
  
  
  for(j in 1:p)
  { ### Middle 90% CI ######
    if(temp_TRUE[j] <= temp_up[j] && temp_TRUE[j] >= temp_low[j])
    {Coverage_Normal_middle[iter_no,j] <- 1}
    
    ##  Upper 90% CI #######
    if(temp_TRUE[j] >= temp_low_10[j])
    {Coverage_Normal_upper[iter_no,j] <- 1}
    
    ##  Lower 90% CI #######
    if(temp_TRUE[j] <= temp_up_90[j])
    {Coverage_Normal_lower[iter_no,j] <- 1}
  }
  Width_Normal[iter_no,] <- as.array(temp_up - temp_low)
  
  
  #########################################################
  ######### Moulton Bootstrap logistic regression #########
  #########################################################
  
  DDD <- bounds_Moulton(X,Y,num_BOOT,x0,iter_no)
  
  temp_low_alt <- DDD$lower_bound_alt
  temp_up_alt <- DDD$upper_bound_alt
  temp_low_10_alt <- DDD$lower_bound_10_alt
  temp_up_90_alt <- DDD$upper_bound_90_alt
  H_star_Moulton_alt_whole_05_10_90_95 <- DDD$H_star_Moulton_alt_whole_05_10_90_95
  
  temp_TRUE <- Beta
  
  H_hat_whole_Moulton_alt <- (sqrt(n)*t(t(Beta_hat - Beta)))
  H_hat_norm <- sum(H_hat_whole_Moulton_alt^2)
  
  ######Coverage of H_n norm ##########################
  
  if(H_hat_norm <= H_star_Moulton_alt_whole_05_10_90_95[4] && H_hat_norm>= H_star_Moulton_alt_whole_05_10_90_95[1])
  {Coverage_Moulton_alt_Hnorm_mid_up_low[iter_no,1] <- 1}
  
  if(H_hat_norm >= H_star_Moulton_alt_whole_05_10_90_95[2])
  {Coverage_Moulton_alt_Hnorm_mid_up_low[iter_no,2] <- 1}
  
  if(H_hat_norm <= H_star_Moulton_alt_whole_05_10_90_95[3])
  {Coverage_Moulton_alt_Hnorm_mid_up_low[iter_no,3] <- 1}
  
  ## Middle 90% CI ######
  for(j in 1:p)
  {if(temp_TRUE[j] <= temp_up_alt[j] && temp_TRUE[j] >= temp_low_alt[j])
  {Coverage_Moulton_alt_middle[iter_no,j] <- 1}
    
    ##  Upper 90% CI #######
    if(temp_TRUE[j] >= temp_low_10_alt[j])
    {Coverage_Moulton_alt_upper[iter_no,j] <- 1}
    
    ##  Lower 90% CI #######
    if(temp_TRUE[j] <= temp_up_90_alt[j])
    {Coverage_Moulton_alt_lower[iter_no,j] <- 1}
  }
  
  Width_Moulton_alt[iter_no,] <- as.array(DDD$width_alt)
  
  #########################################################
  ######### Claeskens Bootstrap logistic regression #########
  #########################################################
  
  DDD <- bounds_Claeskens(X,Y,num_BOOT,x0,iter_no)
  temp_low <- DDD$lower_bound
  temp_up <- DDD$upper_bound
  temp_low_10 <- DDD$lower_bound_10
  temp_up_90 <- DDD$upper_bound_90
  H_3star_whole_05_10_90_95 <- DDD$H_3star_whole_05_10_90_95
  
  
  temp_low_alt <- DDD$lower_bound_alt
  temp_up_alt <- DDD$upper_bound_alt
  temp_low_10_alt <- DDD$lower_bound_10_alt
  temp_up_90_alt <- DDD$upper_bound_90_alt
  H_4star_whole_05_10_90_95 <- DDD$H_4star_whole_05_10_90_95
  
  temp_TRUE <- Beta
  
  H_hat_whole_Claeskens <- (sqrt(n)*t(t(Beta_hat - Beta)))
  H_hat_norm <- sum(H_hat_whole_Claeskens^2)
  
  ######Coverage of H_n norm ##########################
  
  if(H_hat_norm <= H_3star_whole_05_10_90_95[4] && H_hat_norm>= H_3star_whole_05_10_90_95[1])
  {Coverage_Claeskens_Hnorm_mid_up_low[iter_no,1] <- 1}
  
  if(H_hat_norm >= H_3star_whole_05_10_90_95[2])
  {Coverage_Claeskens_Hnorm_mid_up_low[iter_no,2] <- 1}
  
  if(H_hat_norm <= H_3star_whole_05_10_90_95[3])
  {Coverage_Claeskens_Hnorm_mid_up_low[iter_no,3] <- 1}
  
  
  
  if(H_hat_norm <= H_4star_whole_05_10_90_95[4] && H_hat_norm>= H_4star_whole_05_10_90_95[1])
  {Coverage_Claeskens_alt_Hnorm_mid_up_low[iter_no,1] <- 1}
  
  if(H_hat_norm >= H_4star_whole_05_10_90_95[2])
  {Coverage_Claeskens_alt_Hnorm_mid_up_low[iter_no,2] <- 1}
  
  if(H_hat_norm <= H_4star_whole_05_10_90_95[3])
  {Coverage_Claeskens_alt_Hnorm_mid_up_low[iter_no,3] <- 1}
  
  
  ### Middle 90% CI ######
  for(j in 1:p)
  {if(temp_TRUE[j] <= temp_up[j] && temp_TRUE[j] >= temp_low[j])
  {Coverage_Claeskens_middle[iter_no,j] <- 1}
    
    
    ##  Upper 90% CI #######
    if(temp_TRUE[j] >= temp_low_10[j])
    {Coverage_Claeskens_upper[iter_no,j] <- 1}
    
    ##  Lower 90% CI #######
    if(temp_TRUE[j] <= temp_up_90[j])
    {Coverage_Claeskens_lower[iter_no,j] <- 1}
    
    
    
    if(temp_TRUE[j] <= temp_up_alt[j] && temp_TRUE[j] >= temp_low_alt[j])
    {Coverage_Claeskens_alt_middle[iter_no,j] <- 1}
    
    
    ##  Upper 90% CI #######
    if(temp_TRUE[j] >= temp_low_10_alt[j])
    {Coverage_Claeskens_alt_upper[iter_no,j] <- 1}
    
    ##  Lower 90% CI #######
    if(temp_TRUE[j] <= temp_up_90_alt[j])
    {Coverage_Claeskens_alt_lower[iter_no,j] <- 1}
  }
  
  
  Width_Claeskens[iter_no,] <- as.array(DDD$width)
  Width_Claeskens_alt[iter_no,] <- as.array(DDD$width_alt)
  
  
}


### min and max Beta index #######

min_idx <- which(abs(Beta) == min(abs(Beta)))
max_idx <- which(abs(Beta) == max(abs(Beta)))
##################################


mean(Width_Pebble)
mean(Coverage_Pebble_middle)
mean(Coverage_Pebble_upper)
mean(Coverage_Pebble_lower)
colMeans(Coverage_Pebble_Hnorm_mid_up_low)

mean(Width_Normal)
mean(Coverage_Normal_middle)
mean(Coverage_Normal_upper)
mean(Coverage_Normal_lower)
colMeans(Coverage_Normal_Hnorm_mid_up_low)

mean(Width_Moulton_alt)
mean(Coverage_Moulton_alt_middle)
mean(Coverage_Moulton_alt_upper)
mean(Coverage_Moulton_alt_lower)
colMeans(Coverage_Moulton_alt_Hnorm_mid_up_low)

mean(Width_Claeskens)
mean(Coverage_Claeskens_middle)
mean(Coverage_Claeskens_upper)
mean(Coverage_Claeskens_lower)
colMeans(Coverage_Claeskens_Hnorm_mid_up_low)

mean(Width_Claeskens_alt)
mean(Coverage_Claeskens_alt_middle)
mean(Coverage_Claeskens_alt_upper)
mean(Coverage_Claeskens_alt_lower)
colMeans(Coverage_Claeskens_alt_Hnorm_mid_up_low)


if(Num_of_exp == 1)
{Summary_mat <- matrix(0,5,4)

Summary_mat[1,] <- c(mean(Width_Pebble),mean(Coverage_Pebble_middle),mean(Coverage_Pebble_upper),mean(Coverage_Pebble_lower))
Summary_mat[2,] <- c(mean(Width_Normal),mean(Coverage_Normal_middle),mean(Coverage_Normal_upper),mean(Coverage_Normal_lower))
Summary_mat[3,] <- c(mean(Width_Moulton_alt),mean(Coverage_Moulton_alt_middle),mean(Coverage_Moulton_alt_upper),mean(Coverage_Moulton_alt_lower))
Summary_mat[4,] <- c(mean(Width_Claeskens),mean(Coverage_Claeskens_middle),mean(Coverage_Claeskens_upper),mean(Coverage_Claeskens_lower))
Summary_mat[5,] <- c(mean(Width_Claeskens_alt),mean(Coverage_Claeskens_alt_middle),mean(Coverage_Claeskens_alt_upper),
                     mean(Coverage_Claeskens_alt_lower))

rownames(Summary_mat) <- c("PEBBLE","Normal","Moulton","Claeskens 1", "Claeskens 2")
colnames(Summary_mat) <- c("Width","Mid-cvrg.","Up-cvrg.","Low-cvrg.")

Summary_mat_FINAL <- round(Summary_mat,3)

print(Summary_mat_FINAL)
}


###############################################################

if(Num_of_exp > 2)
{Req_array_Pebble <- c(colMeans(Coverage_Pebble_Hnorm_mid_up_low),colMeans(Width_Pebble)[min_idx],
                       sd(Width_Pebble[,min_idx]),
                       colMeans(Coverage_Pebble_middle)[min_idx],colMeans(Coverage_Pebble_upper)[min_idx],
                       colMeans(Coverage_Pebble_lower)[min_idx],colMeans(Width_Pebble)[max_idx],
                       sd(Width_Pebble[,max_idx]),
                       colMeans(Coverage_Pebble_middle)[max_idx],colMeans(Coverage_Pebble_upper)[max_idx],
                       colMeans(Coverage_Pebble_lower)[max_idx],mean(Width_Pebble),
                       sd(Width_Pebble),
                       mean(Coverage_Pebble_middle), mean(Coverage_Pebble_upper),
                       mean(Coverage_Pebble_lower))

Req_array_Normal <- c(colMeans(Coverage_Normal_Hnorm_mid_up_low),colMeans(Width_Normal)[min_idx],
                      sd(Width_Normal[,min_idx]),
                      colMeans(Coverage_Normal_middle)[min_idx],colMeans(Coverage_Normal_upper)[min_idx],
                      colMeans(Coverage_Normal_lower)[min_idx],colMeans(Width_Normal)[max_idx],
                      sd(Width_Normal[,max_idx]),
                      colMeans(Coverage_Normal_middle)[max_idx],colMeans(Coverage_Normal_upper)[max_idx],
                      colMeans(Coverage_Normal_lower)[max_idx],mean(Width_Normal),
                      sd(Width_Normal),
                      mean(Coverage_Normal_middle), mean(Coverage_Normal_upper),
                      mean(Coverage_Normal_lower))

Req_array_Moulton_alt <- c(colMeans(Coverage_Moulton_alt_Hnorm_mid_up_low),colMeans(Width_Moulton_alt)[min_idx],
                           sd(Width_Moulton_alt[,min_idx]),
                           colMeans(Coverage_Moulton_alt_middle)[min_idx],colMeans(Coverage_Moulton_alt_upper)[min_idx],
                           colMeans(Coverage_Moulton_alt_lower)[min_idx],colMeans(Width_Moulton_alt)[max_idx],
                           sd(Width_Moulton_alt[,max_idx]),
                           colMeans(Coverage_Moulton_alt_middle)[max_idx],colMeans(Coverage_Moulton_alt_upper)[max_idx],
                           colMeans(Coverage_Moulton_alt_lower)[max_idx],mean(Width_Moulton_alt),
                           sd(Width_Moulton_alt),
                           mean(Coverage_Moulton_alt_middle), mean(Coverage_Moulton_alt_upper),
                           mean(Coverage_Moulton_alt_lower))

Req_array_Claeskens <- c(colMeans(Coverage_Claeskens_Hnorm_mid_up_low),colMeans(Width_Claeskens)[min_idx],
                         sd(Width_Claeskens[,min_idx]),
                         colMeans(Coverage_Claeskens_middle)[min_idx],colMeans(Coverage_Claeskens_upper)[min_idx],
                         colMeans(Coverage_Claeskens_lower)[min_idx],colMeans(Width_Claeskens)[max_idx],
                         sd(Width_Claeskens[,max_idx]),
                         colMeans(Coverage_Claeskens_middle)[max_idx],colMeans(Coverage_Claeskens_upper)[max_idx],
                         colMeans(Coverage_Claeskens_lower)[max_idx],mean(Width_Claeskens),
                         sd(Width_Claeskens),
                         mean(Coverage_Claeskens_middle), mean(Coverage_Claeskens_upper),
                         mean(Coverage_Claeskens_lower))


Req_array_Claeskens_alt <- c(colMeans(Coverage_Claeskens_alt_Hnorm_mid_up_low),colMeans(Width_Claeskens_alt)[min_idx],
                             sd(Width_Claeskens_alt[,min_idx]),
                             colMeans(Coverage_Claeskens_alt_middle)[min_idx],colMeans(Coverage_Claeskens_alt_upper)[min_idx],
                             colMeans(Coverage_Claeskens_alt_lower)[min_idx],colMeans(Width_Claeskens_alt)[max_idx],
                             sd(Width_Claeskens_alt[,max_idx]),
                             colMeans(Coverage_Claeskens_alt_middle)[max_idx],colMeans(Coverage_Claeskens_alt_upper)[max_idx],
                             colMeans(Coverage_Claeskens_alt_lower)[max_idx],mean(Width_Claeskens_alt),
                             sd(Width_Claeskens_alt),
                             mean(Coverage_Claeskens_alt_middle), mean(Coverage_Claeskens_alt_upper),
                             mean(Coverage_Claeskens_alt_lower))

SUMMARY_MAT <-  matrix(NA,5, length(Req_array_Pebble))

SUMMARY_MAT[1,] <- Req_array_Pebble
SUMMARY_MAT[2,] <- Req_array_Normal
SUMMARY_MAT[3,] <- Req_array_Moulton_alt
SUMMARY_MAT[4,] <- Req_array_Claeskens
SUMMARY_MAT[5,] <- Req_array_Claeskens_alt

rownames(SUMMARY_MAT) <- c("PEBBLE","Normal","Moulton","Claeskens 1", "Claeskens 2")

colnames(SUMMARY_MAT) <- c("H_mid_cvrg","H_up_cvrg","H_low_cvrg","min_id_wide","min_id_wide_sd",
                           "min_id_mid_cvrg", "min_id_up_cvrg","min_id_low_cvrg",
                           "max_id_wide","max_id_wide_sd",
                           "max_id_mid_cvrg", "max_id_up_cvrg","max_id_low_cvrg",
                           "avg_wide","avg_wide_sd",
                           "avg_mid_cvrg", "avg_up_cvrg","avg_low_cvrg")

SUMMARY_MAT_approx <- round(SUMMARY_MAT,3)

print(SUMMARY_MAT_approx[,c(1,4,6,7,8,9,11,12,13)])

write.csv(SUMMARY_MAT_approx, paste0("SummaryALL_mat_nBOOT_",num_BOOT,"_nexp_",Num_of_exp,"_n_",n,"_p_",p,".csv"))
}
# Width_Pebble
# Coverage_Pebble_middle
# Width_Claeskens
# Coverage_Claeskens_middle



# Width_Pebble
# Sigma_hat_diag
# H_star_j_95
# H_star_j_05
# Sigma_hat
# M_hat
# L_hat_inv
# temp_mat_sum
# L_hat
# nearPD(L_hat)


# z_norm_star_record_ALL 
z_norm_ALL 