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
/*
 * Initial computation of Nesterov-Todd scalings.
 * ntsc_n is used for nonlinear constraints.
 * ntsc_l is used for linear cone constraints.
 * ntsc_p is used for second-order cone constraints.
 * ntsc_s is used for positive semidefinite cone constraints.
*/
std::map<std::string,arma::mat> ntsc_n(arma::mat s, arma::mat z){
  std::map<std::string,arma::mat> W;
  int n = s.n_rows;
  arma::mat dnl(n, 1), dnli(n, 1), lambda(n, 1); 
  for(int i = 0; i < n; i++){
    dnl.at(i, 0) = sqrt(s.at(i, 0) / z.at(i, 0));
    dnli.at(i, 0) = sqrt(z.at(i, 0) / s.at(i, 0));
    lambda.at(i, 0) = sqrt(s.at(i, 0) * z.at(i, 0));
  }
  W["dnl"] = dnl;
  W["dnli"] = dnli;
  W["lambda"] = lambda;

  return W;
}
std::map<std::string,arma::mat> ntsc_l(arma::mat s, arma::mat z){
  std::map<std::string,arma::mat> W;
  int n = s.n_rows;
  arma::mat d(n, 1), di(n, 1), lambda(n, 1); 
  for(int i = 0; i < n; i++){
    d.at(i, 0) = sqrt(s.at(i, 0) / z.at(i, 0));
    di.at(i, 0) = sqrt(z.at(i, 0) / s.at(i, 0));
    lambda.at(i, 0) = sqrt(s.at(i, 0) * z.at(i, 0));
  }
  W["d"] = d;
  W["di"] = di;
  W["lambda"] = lambda;

  return W;
}
std::map<std::string,arma::mat> ntsc_p(arma::mat s, arma::mat z){
  std::map<std::string,arma::mat> W;
  int n = s.n_rows;
  arma::mat v(n,1), beta(1,1), lambda(n,1);
  double aa, bb, cc, dd, szdot;

  aa = jnrm2_p(s);
  bb = jnrm2_p(z);
  beta.at(0,0) = sqrt(aa / bb);
  szdot = sdot_nlp(s, z);
  cc = sqrt((szdot / aa / bb + 1.0) / 2.0);
  v = -z / bb;
  v.at(0,0) = -v.at(0,0);
  for(int i = 0; i < n; i++){
    v.at(i,0) += s.at(i,0) / aa;
    v.at(i,0) *= (1.0 / 2.0 / cc);
  }
  v.at(0,0) = v.at(0, 0) + 1.0;
  v *= 1.0 / sqrt(2.0 * v.at(0,0));
  lambda.at(0,0) = cc;
  dd = 2 * cc + s.at(0,0) / aa + z.at(0,0) / bb;
  for(int i = 1; i < n; i++){
    lambda.at(i,0) = s.at(i,0);
    lambda.at(i,0) = (cc + z.at(0,0) / bb) / dd / aa * lambda.at(i,0);
    lambda.at(i,0) = (cc + s.at(0,0) / aa) / dd / bb * z.at(i,0) + lambda.at(i,0); 
    lambda.at(i,0) *= sqrt(aa * bb);
  }
  W["v"] = v;
  W["beta"] = beta;
  W["lambda"] = lambda;

  return W;
}
std::map<std::string,arma::mat> ntsc_s(arma::mat s, arma::mat z, int m){
  std::map<std::string,arma::mat> W;
  arma::mat sc, zc, szc, U, V;
  arma::vec l;
  arma::mat r, rti, lambda;

  s.reshape(m,m);
  z.reshape(m,m);
  arma::chol(sc, s);
  arma::chol(zc, z);
  szc = zc * sc.t();
  svd(U, l, V, szc);
  r = zc.i() * U * diagmat(sqrt(l));
  rti = zc.t() * U * diagmat(1.0 / sqrt(l));
  lambda = diagmat(l);
  lambda.reshape(m * m, 1);
  W["r"] = r;
  W["rti"] = rti;
  W["lambda"] = lambda;

  return W;
}
/*
 * Updating Nesterov-Todd scalings.
 * ntsu_n is used for nonlinear constraints.
 * ntsu_l is used for linear cone constraints.
 * ntsu_p is used for second-order cone constraints.
 * ntsu_s is used for positive semidefinite cone constraints.
*/
std::map<std::string,arma::mat> ntsu_n(std::map<std::string,arma::mat> W, arma::mat s, arma::mat z){
  int n = s.n_rows;
  double ssqrt, zsqrt;

  for(int i = 0; i < n; i++){
    ssqrt = sqrt(s.at(i, 0));
    zsqrt = sqrt(z.at(i, 0));
    W["dnl"].at(i, 0) = W["dnl"].at(i, 0) * ssqrt / zsqrt;
    W["dnli"].at(i, 0) = 1.0 / W["dnl"].at(i, 0);
    W["lambda"].at(i, 0) = ssqrt * zsqrt;
  }

  return W;
}
std::map<std::string,arma::mat> ntsu_l(std::map<std::string,arma::mat> W, arma::mat s, arma::mat z){
  int n = s.n_rows;
  double ssqrt, zsqrt;

  for(int i = 0; i < n; i++){
    ssqrt = sqrt(s.at(i, 0));
    zsqrt = sqrt(z.at(i, 0));
    W["d"].at(i, 0) = W["d"].at(i, 0) * ssqrt / zsqrt;
    W["di"].at(i, 0) = 1.0 / W["d"].at(i, 0);
    W["lambda"].at(i, 0) = ssqrt * zsqrt;
  }

  return W;
}
std::map<std::string,arma::mat> ntsu_p(std::map<std::string,arma::mat> W, arma::mat s, arma::mat z){
  int n = s.n_rows;
  double aa, bb, cc, dd, vs, vz, vq, vu, wk0;

  aa = jnrm2_p(s);
  bb = jnrm2_p(z);
  s /=  aa;
  z /=  bb;
  cc = sqrt((1 + arma::dot(s, z)) / 2.0);
  vs = arma::dot(W["v"], s);
  vz = W["v"].at(0,0) * z.at(0,0);
  for(int i = 1; i < n; i++){
    vz -= W["v"].at(i,0) * z.at(i,0);
  }
  vq = (vs + vz) / 2.0 / cc;
  vu = vs - vz;
  W["lambda"].at(0,0) = cc;
  wk0 = 2 * W["v"].at(0,0) * vq - (s.at(0,0) + z.at(0,0)) / 2.0 / cc;
  dd = (W["v"].at(0,0) * vu - s.at(0,0) / 2.0 + z.at(0,0) / 2.0) / (wk0 + 1.0);
  for(int i = 1; i < n; i++){
    W["lambda"].at(i,0) = W["v"].at(i,0);
    W["lambda"].at(i,0) *= 2.0 * (-dd * vq + 0.5 * vu); 
    W["lambda"].at(i,0) += 0.5 * (1.0 - dd / cc) * s.at(i,0); 
    W["lambda"].at(i,0) += 0.5 * (1.0 + dd / cc) * z.at(i,0); 
    W["lambda"].at(i,0) *= sqrt(aa * bb);
  }
  W["v"] *= 2.0 * vq;
  W["v"].at(0,0) -= s.at(0,0) / 2.0 / cc;
  for(int i = 1; i < n; i++){
    W["v"].at(i,0) += 0.5 / cc * s.at(i,0);
  }
  W["v"] += -0.5 / cc * z;
  W["v"].at(0,0) += 1.0;
  W["v"] *= 1.0 / sqrt(2.0 * W["v"].at(0, 0));
  W["beta"].at(0,0) *= sqrt(aa / bb);

  return W;
}
std::map<std::string,arma::mat> ntsu_s(std::map<std::string,arma::mat> W, arma::mat s, arma::mat z, int m){
  arma::mat zts, U, V, DiagL, lu;
  arma::vec l;

  s.reshape(m,m);
  z.reshape(m,m);
  zts = z.t() * s;
  svd(U, l, V, zts);
  DiagL = diagmat(1.0 / sqrt(l));
  W["r"] = W["r"] * s * V * DiagL;
  W["rti"] = W["rti"] * z * U * DiagL;
  W["lambda"] = diagmat(l);
  W["lambda"].reshape(m * m, 1);

  return W;
}
/*
 * Scaling of vector in S by log-barrier function.
 * sslb_nl is used for nonlinear and linear cone constraints.
 * sslb_p is used for second-order cone constraints.
 * sslb_s is used for positive semidefinite cone constraints.
*/
arma::mat sslb_nl(arma::mat s, arma::mat lambda, bool invers){
  int n = s.n_rows;

  if(invers){
    for(int i = 0; i < n; i++){
      s.at(i,0) *= lambda.at(i,0);
    }
  } else {
    for(int i = 0; i < n; i++){
      s.at(i,0) /= lambda.at(i, 0);
    }
  }

  return s;
}
arma::mat sslb_p(arma::mat s, arma::mat lambda, bool invers){
  int n = s.n_rows;
  double a, cc, lx, s0;

  a = jnrm2_p(lambda);
  if(invers){
    lx = arma::dot(lambda, s) / a;
  } else {
    lx = jdot_p(lambda, s) / a;
  }
  s0 = s.at(0,0);
  s.at(0,0) = lx;
  cc = (lx + s0) / (lambda.at(0,0) / a + 1) / a;
  if(invers == false){
    cc = -cc;
    a = 1 / a;
  }
  for(int i = 1; i < n; i++){
    s.at(i,0) = cc * lambda.at(i,0) + s.at(i,0);
  }
  s = a * s; 

  return s;
}
arma::mat sslb_s(arma::mat s, arma::mat lambda, bool invers, int m){
  arma::vec ld, ls;

  s.reshape(m,m);
  lambda.reshape(m,m);
  ld = lambda.diag();
  for(int i = 0; i < m; i++){
    ls = sqrt(ld(i)) * sqrt(ld);
    if(invers){
      s.col(i) = s.col(i) % ls;
    } else {
      s.col(i) = s.col(i) / ls;
    }
  }
  s.reshape(m * m, 1);

  return s;
}
/*
 * Scaling of vector in S by Nesterov-Todd function.
 * ssnt_n is used for nonlinear constraints.
 * ssnt_l is used for linear cone constraints.
 * ssnt_p is used for second-order cone constraints.
 * ssnt_s is used for positive semidefinite cone constraints.
*/
arma::mat ssnt_n(arma::mat s, std::map<std::string,arma::mat> W, bool invers){
  arma::mat w;
  int m = s.n_rows;
  int n = s.n_cols;

  if(invers){
    w = W["dnli"];
  } else {
    w = W["dnl"];
  }
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      s.at(i, j) *= w.at(i, 0);
    }
  }

  return s;
}
arma::mat ssnt_l(arma::mat s, std::map<std::string,arma::mat> W, bool invers){
  arma::mat w;
  int m = s.n_rows;
  int n = s.n_cols;

  if(invers){
    w = W["di"];
  } else {
    w = W["d"];
  }
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      s.at(i, j) *= w.at(i, 0);
    }
  }

  return s;
}
arma::mat ssnt_p(arma::mat s, std::map<std::string,arma::mat> W, bool invers){
  double a;
  arma::mat w;

  if(invers){
    s.row(0) = -s.row(0);
  }
  w = s.t() * W["v"];
  s.row(0) = -s.row(0);
  s = 2 * W["v"] * w.t() + s;
  if(invers){
    s.row(0) = -s.row(0);
    a = 1 / W["beta"].at(0,0);
  } else {
    a = W["beta"].at(0,0);
  }
  s = a * s;

  return s;
}
arma::mat ssnt_s(arma::mat s, std::map<std::string,arma::mat> W, bool invers, bool transp){
  int n, m;
  bool tt;
  arma::mat w, S, a, ans;

  if(invers){
    w = W["rti"];
    tt = transp;
  } else {
    w = W["r"];
    if(transp){
      tt = false;
    } else {
      tt = true;
    }
  }
  m = w.n_cols;
  n = s.n_cols;
  for(int i = 0; i < n; i++){
    S = s.col(i);
    S.reshape(m,m);
    S.diag() = 0.5 * S.diag();
    for(int k = 0; k < m; k++){
      for(int r = 0; r < k; r++){
	S.at(r,k) = 0.0;
      }
    } 
    if(tt){
      a = S * w;
      ans = w.t() * a + a.t() * w;
    } else {
      a = w * S;
      ans = w * a.t() + a * w.t();
    }
    ans.reshape(m * m, 1);
    s.col(i) = ans;
  }

  return s;
}

