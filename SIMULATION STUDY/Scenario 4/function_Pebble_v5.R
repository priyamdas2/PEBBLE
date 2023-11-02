########### Defining Regular logistic likelihood function ONLY beta ###########
negloglikelogistic_ORIGINAL_BETA <- function(Beta){
  return(-loglikelogistic_ORIGINAL(Beta, Y, X))
}


########### Defining Bootstrap logistic likelihood function ONLY beta ###########
negloglikelogistic_BETA <- function(Beta){
  return(-loglikelogistic(Beta, Beta_hat, Y, X, G, mu_G))
}

###############################################################################
########### Defining T_star  ##################################################
###############################################################################
T_star <- function(X,Y,X_0,Beta_hat,Beta_star, z_norm_star,G,mu_G, b)
{n <- dim(X)[1]
p <- dim(X)[2]

temp_prod <- exp(X%*%Beta_star)
coef_i <- temp_prod/((1+temp_prod)^2)     # For hat(L_n)^*

temp_prod_2 <- exp(X%*%Beta_hat)
p_hat <- temp_prod_2/(1+temp_prod_2)      # For hat(L_n)^*

temp_mat_sum <- matrix(0, p, p)
temp_mat_sum_2 <- matrix(0, p, p)
for(i in 1:n)
{prod_temp <- t(t(X[i,]))%*%X[i,]
temp_mat_sum <- temp_mat_sum + prod_temp*coef_i[i]   # For hat(L_n)^*
temp_mat_sum_2 <- temp_mat_sum_2 + 
  prod_temp*((Y[i] - p_hat[i])^2)*((G[i] - mu_G)^2)/(mu_G^2)}  # For hat(M_n)^*

L_star <- temp_mat_sum/n
# L_star <- temp_mat_sum/n + diag(p)/sqrt(n)
L_star_obj <- nearPD(L_star)
L_star <- L_star_obj$mat
L_star_inv <- solve(L_star)

M_star <- temp_mat_sum_2/n

Sigma_star <- L_star_inv%*%M_star%*%L_star_inv
Sigma_star_obj <- nearPD(Sigma_star)
Sigma_star <- Sigma_star_obj$mat

part_1 <- exp(-sum(X_0*Beta_star))*(t(X_0)%*%Sigma_star%*%X_0)^(-0.5)
part_2 <- sqrt(n)*(exp(sum(X_0*Beta_star)) - exp(sum(X_0*Beta_hat))) +
          b*exp(sum(X_0*Beta_star))*t(X_0)%*%L_star_inv%*%t(t(z_norm_star))
T_star <- part_1*part_2
return(T_star)
}  
  



###############################################################################
########### Defining H_star_compwise ##################################################
###############################################################################
H_star <- function(X,Y,Beta_hat,Beta_star, z_norm_star,G,mu_G, b)
{n <- dim(X)[1]
p <- dim(X)[2]

temp_prod <- exp(X%*%Beta_star)
coef_i <- temp_prod/((1+temp_prod)^2)     # For hat(L_n)^*

temp_prod_2 <- exp(X%*%Beta_hat)
p_hat <- temp_prod_2/(1+temp_prod_2)      # For hat(L_n)^*

temp_mat_sum <- matrix(0, p, p)
temp_mat_sum_2 <- matrix(0, p, p)
for(i in 1:n)
{prod_temp <- t(t(X[i,]))%*%X[i,]
temp_mat_sum <- temp_mat_sum + prod_temp*coef_i[i]   # For hat(L_n)^*
temp_mat_sum_2 <- temp_mat_sum_2 + 
  prod_temp*((Y[i] - p_hat[i])^2)*((G[i] - mu_G)^2)/(mu_G^2)}  # For hat(M_n)^*

L_star <- temp_mat_sum/n
# L_star <- temp_mat_sum/n + diag(p)/sqrt(n)
L_star_obj <- nearPD(L_star)
L_star <- L_star_obj$mat
L_star_inv <- solve(L_star)

M_star <- temp_mat_sum_2/n


# \Sigma^*_j = (j,j) - th element of (L^*)^{-1} * (M^*) * (L^*)^{-1}

Sigma_star <- L_star_inv%*%M_star%*%L_star_inv
Sigma_star_obj <- nearPD(Sigma_star)
Sigma_star <- Sigma_star_obj$mat
Sigma_star_inv <- solve(Sigma_star)
Sigma_star_pow_minus_half <- sqrtm(Sigma_star_inv)

H_star_whole <- Sigma_star_pow_minus_half%*%(sqrt(n)*t(t(Beta_star - Beta_hat))+b*L_star_inv%*%t(t(z_norm_star)))
H_star_norm <- sum(H_star_whole^2)

Sig_diag_pow_minus_half <- array(NA,p)
H_star_compwise <- matrix(NA,p,1)

for(jj in 1:p)
{Sig_diag_pow_minus_half[jj] <- (Sigma_star[jj,jj])^(-1/2)
last_part <- b*L_star_inv[jj,]%*%t(t(z_norm_star))
H_star_compwise[jj,1] <- Sig_diag_pow_minus_half[jj]*(sqrt(n)*(Beta_star[jj] - Beta_hat[jj]) + last_part)}

return(list(H_star_compwise = H_star_compwise, H_star_norm = H_star_norm))
}
###############################################################################
###############################################################################
###############################################################################




####################################################################################
###############  Define bounds_Pebble ################################################
####################################################################################

bounds_Pebble <- function(X,Y,X_0,b, num_BOOT, beta_a, beta_b,x0,iter_no,Sigma_hat_diag,Norm_var)
{n <- dim(X)[1]
p <- dim(X)[2]
bound_max <- 10
#### Finding Beta_hat and L_hat ##############
negloglikelogistic_ORIGINAL_BETA <- function(Beta){
  return(-loglikelogistic_ORIGINAL(Beta, Y, X))
}
x0 <- array(1,p)

dd <- optimr(x0,negloglikelogistic_ORIGINAL_BETA)
Beta_hat <- dd$par
#####################################################

H_star_MATRIX <- matrix(NA,num_BOOT,p)
H_star_norms <- array(NA,num_BOOT)
z_norm_star_record <- matrix(NA,num_BOOT,p)
T_star_values <- array(NA,num_BOOT)

diff_in_beta_approx_star_MATRIX <- matrix(NA,num_BOOT,p)

Sigma_star_pow_neg_half_ARRAY <- array(NA,num_BOOT)
Sigma_star_MAT <- matrix(0,num_BOOT,p)


for(i in 1:num_BOOT)
{ #print(i)
  #message(paste('PEBBLE Bootstrap iteration number =',i,',iter_no = ',iter_no))
  #### Finding Beta_star    #########
  G <- rbeta(n, beta_a, beta_b)
  mu_G <- beta_a/(beta_a+beta_b)
  
  ########### Defining Regular logistic likelihood function ONLY beta ###########
  negloglikelogistic_BETA <- function(Beta){
    return(-loglikelogistic(Beta, Beta_hat, Y, X, G, mu_G))
  }
  ###############################################################################
  dd <- optimr(x0,negloglikelogistic_BETA)
  Beta_star <- dd$par 
  
  z_norm_star <- rnorm(p,0,Norm_var)
  z_norm_star_record[i,] <- z_norm_star
  
  ##### Finding H,ub,lb, etc #########
  DD <- H_star(X,Y,Beta_hat,Beta_star, z_norm_star,G,mu_G, b)
  H_star_MATRIX[i,] <- as.array(DD$H_star_compwise)
  H_star_norms[i] <- DD$H_star_norm
  #diff_in_beta_approx_star_MATRIX[i,] <- (H_star_MATRIX[i,]*sqrt(Sigma_hat_diag))/sqrt(n)
  
  ###### Finding T quantiles #########
  EE <- T_star(X,Y,X_0,Beta_hat,Beta_star, z_norm_star,G,mu_G, b)
  T_star_values[i] <- EE
}



row_maxs <- apply(abs(H_star_MATRIX), 1, max)
to_be_selected <- which(row_maxs <= bound_max)

H_star_j_05 <- array(NA,p)
H_star_j_10 <- array(NA,p)
H_star_j_90 <- array(NA,p)
H_star_j_95 <- array(NA,p)

H_star_whole_05_10_90_95 <- array(NA,4)

if(length(to_be_selected) >= 20)
{index_q_05 <- array(NA,p)
index_q_10 <- array(NA,p)
index_q_90 <- array(NA,p)
index_q_95 <- array(NA,p)

H_star_MATRIX_SELECTED <- H_star_MATRIX[to_be_selected,]
H_star_norms_SELECTED <- H_star_norms[to_be_selected]
T_star_values_SELECTED <- T_star_values[to_be_selected]

for(j in 1:p)
{index_q_05[j] <- which.min(abs(H_star_MATRIX_SELECTED[,j] - quantile(H_star_MATRIX_SELECTED[,j], .05)))
index_q_10[j] <- which.min(abs(H_star_MATRIX_SELECTED[,j] - quantile(H_star_MATRIX_SELECTED[,j], .1)))
index_q_90[j] <- which.min(abs(H_star_MATRIX_SELECTED[,j] - quantile(H_star_MATRIX_SELECTED[,j], .90)))
index_q_95[j] <- which.min(abs(H_star_MATRIX_SELECTED[,j] - quantile(H_star_MATRIX_SELECTED[,j], .95)))

H_star_j_05[j] <- H_star_MATRIX_SELECTED[index_q_05[j],j]
H_star_j_10[j] <- H_star_MATRIX_SELECTED[index_q_10[j],j]
H_star_j_90[j] <- H_star_MATRIX_SELECTED[index_q_90[j],j]
H_star_j_95[j] <- H_star_MATRIX_SELECTED[index_q_95[j],j]
}

index_whole <- which.min(abs(H_star_norms_SELECTED - quantile(H_star_norms_SELECTED, .05)))
index_whole[2] <- which.min(abs(H_star_norms_SELECTED - quantile(H_star_norms_SELECTED, .1)))
index_whole[3] <- which.min(abs(H_star_norms_SELECTED - quantile(H_star_norms_SELECTED, .9)))
index_whole[4] <- which.min(abs(H_star_norms_SELECTED - quantile(H_star_norms_SELECTED, .95)))

H_star_whole_05_10_90_95 <- H_star_norms_SELECTED[index_whole]

T_star_quantiles_05_10_90_95 <- c(quantile(T_star_values_SELECTED,0.05), quantile(T_star_values_SELECTED,0.1),
                                  quantile(T_star_values_SELECTED,0.9),quantile(T_star_values_SELECTED,0.95))
}

return(list(H_star_j_05 = H_star_j_05, H_star_j_10 = H_star_j_10, H_star_j_90 = H_star_j_90,
            H_star_j_95 = H_star_j_95, H_star_whole_05_10_90_95 = H_star_whole_05_10_90_95,
            z_norm_star_record = z_norm_star_record,
            T_star_quantiles_05_10_90_95 = T_star_quantiles_05_10_90_95))
}

#########################################################################################
#########################################################################################
#########################################################################################

