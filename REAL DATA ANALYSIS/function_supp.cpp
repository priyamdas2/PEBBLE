#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>


// [[Rcpp::depends(RcppArmadillo)]]

#define crossprod(x) x.t() * x
#define tcrossprod(x) x * x.t()


using namespace arma;

// [[Rcpp::export]]
double ceil_fun(double number) {
  return ceil(number);
}

// [[Rcpp::export]]
double loglikelogistic_ORIGINAL(vec Beta,vec Y, mat X) {
  int n = Y.size();
  int p = X.n_cols;
  vec prod_vec(n);
  vec P_x(n);
  vec loglike(n);
  double total_loglike = 0;


  for(int i = 0; i < n; ++i) {
    prod_vec[i] = 0;
    for(int j = 0; j < p; ++j)
    {prod_vec[i] = prod_vec[i] + X(i,j)*Beta[j];
    }
    P_x[i] = exp(prod_vec[i])/(1+ exp(prod_vec[i]));

    loglike[i] = Y[i]*log(P_x[i]) + (1-Y[i])*log(1-P_x[i]);
    total_loglike = total_loglike + loglike[i];
  }

  return total_loglike;
}

// [[Rcpp::export]]
double loglikelogistic(vec Beta, vec Beta_hat, vec Y, mat X, vec G, double mu_G) {
  int n = Y.size();
  int p = X.n_cols;
  vec prod_vec(n);
  vec prod_vec_hat(n);
  vec P_x(n);
  vec P_x_hat(n);
  vec loglike(n);
  double total_loglike = 0;


  for(int i = 0; i < n; ++i) {
    prod_vec[i] = 0;
    prod_vec_hat[i] = 0;
    for(int j = 0; j < p; ++j)
    {prod_vec[i] = prod_vec[i] + X(i,j)*Beta[j];
      prod_vec_hat[i] = prod_vec_hat[i] + X(i,j)*Beta_hat[j];
    }
    P_x_hat[i] = exp(prod_vec_hat[i])/(1+ exp(prod_vec_hat[i]));

    loglike[i] = (Y[i]-P_x_hat[i])*prod_vec[i]*(G[i] - mu_G) + mu_G*(P_x_hat[i]*prod_vec[i] - log(1+exp(prod_vec[i])));
    total_loglike = total_loglike + loglike[i];
  }

  return total_loglike;
}


// [[Rcpp::export]]
double loglikelogistic2(vec Beta, vec Beta_hat, vec Y, mat X, vec G, double mu_G) {
  int n = Y.size();
  int p = X.n_cols;
  vec prod_vec(n);
  vec prod_vec_hat(n);
  vec P_x(n);
  vec P_x_hat(n);
  vec Semi_P_x_hat(n);
  vec loglike(n);
  double total_loglike = 0;
  
  
  for(int i = 0; i < n; ++i) {
    prod_vec[i] = 0;
    prod_vec_hat[i] = 0;
    for(int j = 0; j < p; ++j)
    {prod_vec[i] = prod_vec[i] + X(i,j)*Beta[j];
      prod_vec_hat[i] = prod_vec_hat[i] + X(i,j)*Beta_hat[j];
    }
    P_x_hat[i] = exp(prod_vec_hat[i])/(1+ exp(prod_vec_hat[i]));
    Semi_P_x_hat[i] = sqrt(exp(prod_vec_hat[i]))/(1+ exp(prod_vec_hat[i]));
    
    loglike[i] = Semi_P_x_hat[i]*prod_vec[i]*(G[i] - mu_G) + mu_G*(P_x_hat[i]*prod_vec[i] - log(1+exp(prod_vec[i])));
    total_loglike = total_loglike + loglike[i];
  }
  
  return total_loglike;
}
