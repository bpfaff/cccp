#include "cccp3.h"

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif
/*
 * Inner product of two vectors in S.
 * snrm2_nls is used for nonlinear, linear and second-order cone constraints.
 * snrm2_p is used for positive semi-definite constraints.
*/
double snrm2_nls(arma::mat s) {
  return norm(s);
}
double snrm2_p(arma::mat s, int m) {
  double ans = 0.0;
  ans = sqrt(sdot_p(s, s, m));
  return ans;
}
