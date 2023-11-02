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
mat H_star(mat X, vec Beta_hat, vec Beta_star, vec z_norm, double c) {
  int n = X.n_rows;
  int p = X.n_cols;
  int p_1 = 4;
  double b;
  double temp;
  mat L_star(p,p);
  mat H_star_final(p,p);
  mat Sigma_star(p,p);
  mat Sigma_star_pow_neg_half(p,p);
  mat temp_diag(p,p);
  mat temp_mat_sum;
  vec coef_i;
  double temp_const;

  mat U;
  vec s;
  mat V;

  if(p_1 > 4)
  {p_1 = p + 1;
   }
  temp = pow(n,-(1/p_1));
  b = pow(c*temp,0.5);

  for(int i = 0; i < p; ++i)
  {temp_const = sum(dot(X.row(i),Beta_star));  // WRONG here, should be X.row(i)
   coef_i[i] = exp(temp_const)/pow((1+exp(temp_const)),2);
   temp_mat_sum = temp_mat_sum + coef_i[i]*(X.row(i).t()*X.row(i));
  }
  L_star = temp_mat_sum/n;
  Sigma_star = inv(L_star);

  // Finding Sigma_star_pow_neg_half
  svd_econ(U, s, V, L_star);
  for(int i = 0; i < p; ++i)
  {temp_diag(i,i) = sqrt(s[i]);
    }

  Sigma_star_pow_neg_half = U*temp_diag;

  H_star_final = Sigma_star_pow_neg_half*(sqrt(n)*(Beta_star.as_col() - Beta_hat.as_col()) + b*Sigma_star*z_norm.as_col());

  return H_star_final;
}




