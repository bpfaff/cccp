#include "cccp3.h"

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Inner product of two vectors in S.
 * sdot_nls is used for nonlinear, linear and second-order cone constraints.
 * sdot_p is used for positive semi-definite constraints.
*/
double sdot_nls(arma::mat s, arma::mat z) {
  double ans = arma::dot(s, z);
  return ans;
}
double sdot_p(arma::mat s, arma::mat z, int m) {
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
 * Inner product of two vectors in S.
 * snrm2_nls is used for nonlinear, linear and second-order cone constraints.
 * snrm2_p is used for positive semi-definite constraints.
*/
double snrm2_nls(arma::mat s) {
  return arma::norm(s);
}
double snrm2_p(arma::mat s, int m) {
  double ans = 0.0;
  ans = sqrt(sdot_p(s, s, m));
  return ans;
}
