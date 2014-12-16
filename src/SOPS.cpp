#include "cccp3.h"

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Inner product of two vectors in S.
 * sdot_nlp is used for nonlinear, linear and second-order cone constraints.
 * sdot_s is used for positive semi-definite constraints.
*/
double sdot_nlp(arma::mat s, arma::mat z) {
  double ans = arma::dot(s, z);
  return ans;
}
double sdot_s(arma::mat s, arma::mat z, int m) {
  double ans = 0.0;
  int n = s.n_rows;

  // squaring and summing diagonal elements
  for(int i = 0; i < n; i += m + 1){
    ans += s.at(i,0) * z.at(i,0);
  }
  // product of lower-diagonal elements, multiplied by two
  for(int i = 0; i < m; i++){
    for(int j = 0; j < (m - 1); j++){
      if(j < i){
	ans += 2.0 * s.at(i + m * j,0) * z.at(i + m * j, 0);
      }
    }
  }
  return ans;
}
/*
 * J-dot product between two vectors related to second-order cone constraints.
 * Returns x' * J * y, whereby J = [1, 0; 0, -I]
 * jdot_p
*/
double jdot_p(arma::mat s, arma::mat z){
  int n = s.n_rows;
  double ans;

  ans = s.at(0,0) * z.at(0,0);
  for(int i = 1; i < n; i++){
    ans -= s(i, 0) * z(i, 0);
  }

  return ans;
}
/*
 * Norm of a vectors in S.
 * snrm2_nlp is used for nonlinear, linear and second-order cone constraints.
 * snrm2_s is used for positive semi-definite constraints.
*/
double snrm2_nlp(arma::mat s) {
  return arma::norm(s);
}
double snrm2_s(arma::mat s, int m) {
  double ans = 0.0;
  ans = sqrt(sdot_s(s, s, m));
  return ans;
}
/*
 * J-norm of a vector related to a second-order cone constraint.
 * Returns the square-root of x' * J * y, whereby J = [1, 0; 0, -I]
 * jnrm2_p
*/
double jnrm2_p(arma::mat s){
  double ans = sqrt(jdot_p(s, s));
  return ans;
}
/*
 * Product between two vectors in S.
 * sprd_nl is used for nonlinear and linear cone constraints.
 * sprd_p is used for second-order cone constraints.
 * sprd_s is used for positive semidefinite cone constraints.
*/
arma::mat sprd_nl(arma::mat s, arma::mat z){
  return s % z;
}
arma::mat sprd_p(arma::mat s, arma::mat z){
  int n = s.n_rows;
  arma::mat ans(n, 1);

  ans(0, 0) = arma::dot(s, z);
  for(int i = 1; i < n; i++){
    ans(i, 0) = s(0, 0) * z(i, 0) + z(0, 0) * s(i, 0);
  }

  return ans;
}
arma::mat sprd_s(arma::mat s, arma::mat z, int m){
  arma::mat ans(m,m);
  s.reshape(m,m);
  z.reshape(m,m);

  ans = s * z;
  ans.reshape(m * m, 1);

  return ans;
}
/*
 * One-element (neutral) with respect to a vector in S.
 * sone_nl is used for nonlinear and linear cone constraints.
 * sone_p is used for second-order cone constraints.
 * sone_s is used for positive semidefinite cone constraints.
*/
arma::mat sone_nl(arma::mat s){
  arma::mat ans(s.n_rows, 1);
  ans.ones();
  return ans;
}
arma::mat sone_p(arma::mat s){
  arma::mat ans(s.n_rows, 1);
  ans.zeros();
  ans.at(0,0) = 1.0;
  return ans;
}
arma::mat sone_s(int m){
  arma::mat ans = arma::eye(m, m);
  ans.reshape(m * m, 1);
  return ans;
}
/*
 * Inverse of product between two vectors in S.
 * sinv_nl is used for nonlinear and linear cone constraints.
 * sinv_p is used for second-order cone constraints.
 * sinv_s is used for positive semidefinite cone constraints.
*/
arma::mat sinv_nl(arma::mat s, arma::mat z){
  arma::mat ans(s.n_rows, 1);
  ans = s / z;
  return ans;
}
arma::mat sinv_p(arma::mat s, arma::mat z){
  int n = s.n_rows;
  arma::mat ans(n,1);
  double aa = jdot_p(z, z);
  double cc = s(0, 0);
  double dd = 0.0;

  for(int i = 1; i < n; i++){
    dd += s(i, 0) * z(i, 0);
  }
  ans(0, 0) = cc * z(0, 0) - dd;
  for(int i = 1; i < n; i++){
    ans(i, 0) = aa / z(0, 0) * s(i, 0);
    ans(i, 0) += (dd / z(0, 0) - cc) * z(i, 0);
  }
  for(int i = 0; i < n; i++){
    ans(i, 0) /= aa;
  }

  return ans;
}
arma::mat sinv_s(arma::mat s, arma::mat z, int m){
  arma::mat ans(m,m);

  s.reshape(m,m);
  z.reshape(m,m);
  ans = inv(z) * s;
  ans.reshape(m * m, 1);

  return ans;
}
/*
 * Determining maximum step-size of a vector in S.
 * smss_nl is used for nonlinear and linear cone constraints.
 * smss_p is used for second-order cone constraints.
 * smss_s is used for positive semidefinite cone constraints.
*/
double smss_nl(arma::mat s){
  return s.min();
}
double smss_p(arma::mat s){
  double ans = 0.0;
  int n = s.n_rows;

  for(int i = 1; i < n; i++){
    ans += s(i, 0) * s(i, 0);
  }
  ans = sqrt(ans);
  ans = ans - s(0, 0);

  return ans;
}
double smss_s(arma::mat s, int m){
  arma::vec eval;
  arma::mat evec;

  s.reshape(m,m);
  arma::eig_sym(eval, evec, s);

  return -eval[0];
}
/*
 * Applying maximum step-size to a vector in S (initial).
 * sams1_nl is used for nonlinear and linear cone constraints.
 * sams1_p is used for second-order cone constraints.
 * sams1_s is used for positive semidefinite cone constraints.
*/
arma::mat sams1_nl(arma::mat s, double alpha){
  arma::mat adj(1,1);

  adj.at(0,0) = 1 + alpha;
  s.each_row() += adj;

  return s; 
}
arma::mat sams1_p(arma::mat s, double alpha){
  s.at(0,0) += 1.0 + alpha;

  return s;
}
arma::mat sams1_s(arma::mat s, double alpha, int m){
  s.reshape(m,m);
  s.diag() = s.diag() + (1 + alpha);
  s.reshape(m * m, 1);

  return s;
}
/*
 * Applying maximum step-size to a vector in S (during iterations).
 * sams1_nl is used for nonlinear and linear cone constraints.
 * sams1_p is used for second-order cone constraints.
 * sams1_s is used for positive semidefinite cone constraints.
*/
arma::mat sams2_nl(arma::mat s, double alpha){
  int n = s.n_rows;

  for(int i = 0; i < n; i++){
    s.at(i, 0) = 1.0 + alpha * s.at(i, 0);
  }

  return s; 
}
arma::mat sams2_p(arma::mat s, double alpha){
  int n = s.n_rows;

  for(int i = 0; i < n; i++){
    s.at(i, 0) = alpha * s.at(i, 0);
  }
  s.at(0, 0) += 1.0;

  return s;
}
arma::mat sams2_s(arma::mat s, double alpha, arma::mat lambda, arma::vec sigma, int m){
  s.reshape(m,m);
  lambda.reshape(m,m);

  for(int i = 0; i < m; i++){
    sigma.at(i) = 1 + alpha * sigma.at(i);
    s.col(i) = s.col(i) * sqrt(sigma.at(i) / lambda.at(i, i)); 
  }
  s.reshape(m * m, 1);

  return s;
}
