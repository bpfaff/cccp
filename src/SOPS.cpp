#include "cccp.h"

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif
using namespace arma;
/*
 * Inner product of two vectors in S.
 * sdot_nlp is used for nonlinear, linear and second-order cone constraints.
 * sdot_s is used for positive semi-definite constraints.
*/
double sdot_nlp(mat s, mat z) {
  double ans = dot(s, z);
  return ans;
}
double sdot_s(mat s, mat z, int m) {
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
double jdot_p(mat s, mat z){
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
double snrm2_nlp(mat s) {
  return norm(s);
}
double snrm2_s(mat s, int m) {
  double ans = 0.0;
  ans = sqrt(sdot_s(s, s, m));
  return ans;
}
/*
 * J-norm of a vector related to a second-order cone constraint.
 * Returns the square-root of x' * J * y, whereby J = [1, 0; 0, -I]
 * jnrm2_p
*/
double jnrm2_p(mat s){
  double ans = sqrt(jdot_p(s, s));
  return ans;
}
/*
 * Product between two vectors in S.
 * sprd_nl is used for nonlinear and linear cone constraints.
 * sprd_p is used for second-order cone constraints.
 * sprd_s is used for positive semidefinite cone constraints.
*/
mat sprd_nl(mat s, mat z){
  return s % z;
}
mat sprd_p(mat s, mat z){
  int n = s.n_rows;
  mat ans(n, 1);

  ans(0, 0) = dot(s, z);
  for(int i = 1; i < n; i++){
    ans(i, 0) = s(0, 0) * z(i, 0) + z(0, 0) * s(i, 0);
  }

  return ans;
}
mat sprd_s(mat s, mat z, int m){
  mat ans(m,m);
  s.reshape(m,m);
  z.reshape(m,m);

  ans = 0.5 * (s * z + z * s);
  ans.reshape(m * m, 1);

  return ans;
}
/*
 * One-element (neutral) with respect to a vector in S.
 * sone_nl is used for nonlinear and linear cone constraints.
 * sone_p is used for second-order cone constraints.
 * sone_s is used for positive semidefinite cone constraints.
*/
mat sone_nl(int m){
  mat ans(m, 1);
  ans.ones();
  return ans;
}
mat sone_p(int m){
  mat ans(m, 1);
  ans.zeros();
  ans.at(0,0) = 1.0;
  return ans;
}
mat sone_s(int m){
  mat ans = eye(m, m);
  ans.reshape(m * m, 1);
  return ans;
}
/*
 * Inverse of product between two vectors in S.
 * sinv_nl is used for nonlinear and linear cone constraints.
 * sinv_p is used for second-order cone constraints.
 * sinv_s is used for positive semidefinite cone constraints.
*/
mat sinv_nl(mat s, mat z){
  mat ans(s.n_rows, 1);
  ans = s / z;
  return ans;
}
mat sinv_p(mat s, mat z){
  int n = s.n_rows;
  mat ans(n,1);
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
mat sinv_s(mat s, mat z, int m){
  mat ans(m,m);

  s.reshape(m,m);
  z.reshape(m,m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      ans.at(i,j) = s.at(i,j) * 2.0 / (z.at(i,i) + z.at(j,j));
	}
  }
  ans.reshape(m * m, 1);

  return ans;
}
/*
 * Determining maximum step-size of a vector in S.
 * smss_nl is used for nonlinear and linear cone constraints.
 * smss_p is used for second-order cone constraints.
 * smss_s is used for positive semidefinite cone constraints.
*/
double smss_nl(mat s){
  return -s.min();
}
double smss_p(mat s){
  double ans = 0.0;
  int n = s.n_rows;

  for(int i = 1; i < n; i++){
    ans += s(i, 0) * s(i, 0);
  }
  ans = sqrt(ans);
  ans = ans - s(0, 0);

  return ans;
}
double smss_s(mat s, int m){
  vec eval;
  mat evec;

  s.reshape(m,m);
  eig_sym(eval, evec, s);

  return -eval.min();
}
/*
 * Applying maximum step-size to a vector in S (initial).
 * sams1_nl is used for nonlinear and linear cone constraints.
 * sams1_p is used for second-order cone constraints.
 * sams1_s is used for positive semidefinite cone constraints.
*/
mat sams1_nl(mat s, double alpha){
  mat adj(1,1);

  adj.at(0,0) = 1 + alpha;
  s.each_row() += adj;

  return s; 
}
mat sams1_p(mat s, double alpha){
  s.at(0,0) += 1.0 + alpha;

  return s;
}
mat sams1_s(mat s, double alpha, int m){
  s.reshape(m,m);
  s.diag() = s.diag() + (1 + alpha);
  s.reshape(m * m, 1);

  return s;
}
/*
 * Applying maximum step-size to a vector in S (during iterations).
 * sams2_nl is used for nonlinear and linear cone constraints.
 * sams2_p is used for second-order cone constraints.
 * sams2_s is used for positive semidefinite cone constraints.
*/
mat sams2_nl(mat s, double alpha){
  int n = s.n_rows;

  for(int i = 0; i < n; i++){
    s.at(i, 0) = 1.0 + alpha * s.at(i, 0);
  }

  return s; 
}
mat sams2_p(mat s, double alpha){
  int n = s.n_rows;

  for(int i = 0; i < n; i++){
    s.at(i, 0) = alpha * s.at(i, 0);
  }
  s.at(0, 0) += 1.0;

  return s;
}
mat sams2_s(mat s, double alpha, mat lambda, vec sigma, int m){
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
std::map<std::string,mat> ntsc_n(mat s, mat z){
  std::map<std::string,mat> W;
  int n = s.n_rows;
  mat dnl(n, 1), dnli(n, 1), lambda(n, 1); 
  for(int i = 0; i < n; i++){
    dnl.at(i, 0) = sqrt(s.at(i, 0) / z.at(i, 0));
    dnli.at(i, 0) = sqrt(z.at(i, 0) / s.at(i, 0));
    lambda.at(i, 0) = sqrt(s.at(i, 0) * z.at(i, 0));
  }
  W.insert(std::pair<std::string,mat>("dnl", dnl));
  W.insert(std::pair<std::string,mat>("dnli", dnli));
  W.insert(std::pair<std::string,mat>("lambda", lambda));

  return W;
}
std::map<std::string,mat> ntsc_l(mat s, mat z){
  std::map<std::string,mat> W;
  int n = s.n_rows;
  mat d(n, 1), di(n, 1), lambda(n, 1); 
  for(int i = 0; i < n; i++){
    d.at(i, 0) = sqrt(s.at(i, 0) / z.at(i, 0));
    di.at(i, 0) = sqrt(z.at(i, 0) / s.at(i, 0));
    lambda.at(i, 0) = sqrt(s.at(i, 0) * z.at(i, 0));
  }
  W.insert(std::pair<std::string,mat>("d", d));
  W.insert(std::pair<std::string,mat>("di", di));
  W.insert(std::pair<std::string,mat>("lambda", lambda));

  return W;
}
std::map<std::string,mat> ntsc_p(mat s, mat z){
  std::map<std::string,mat> W;
  int n = s.n_rows;
  mat v(n,1), beta(1,1), lambda(n,1);
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
  }
  lambda = sqrt(aa * bb) * lambda;
  W.insert(std::pair<std::string,mat>("v", v));
  W.insert(std::pair<std::string,mat>("beta", beta));
  W.insert(std::pair<std::string,mat>("lambda", lambda));

  return W;
}
std::map<std::string,mat> ntsc_s(mat s, mat z, int m){
  std::map<std::string,mat> W;
  mat sc, zc, szc, U, V;
  vec l;
  mat r, rti, lambda;

  s.reshape(m,m);
  z.reshape(m,m);
  chol(sc, s);
  chol(zc, z);
  szc = zc * sc.t();
  svd(U, l, V, szc);
  r = zc.i() * U * diagmat(sqrt(l));
  rti = zc.t() * U * diagmat(1.0 / sqrt(l));
  lambda = diagmat(l);
  lambda.reshape(m * m, 1);
  W.insert(std::pair<std::string,mat>("r", r));
  W.insert(std::pair<std::string,mat>("rti", rti));
  W.insert(std::pair<std::string,mat>("lambda", lambda));

  return W;
}
/*
 * Updating Nesterov-Todd scalings.
 * ntsu_n is used for nonlinear constraints.
 * ntsu_l is used for linear cone constraints.
 * ntsu_p is used for second-order cone constraints.
 * ntsu_s is used for positive semidefinite cone constraints.
*/
std::map<std::string,mat> ntsu_n(std::map<std::string,mat> W, mat s, mat z){
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
std::map<std::string,mat> ntsu_l(std::map<std::string,mat> W, mat s, mat z){
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
std::map<std::string,mat> ntsu_p(std::map<std::string,mat> W, mat s, mat z){
  int n = s.n_rows;
  double aa, bb, cc, dd, vs, vz, vq, vu, wk0;

  aa = jnrm2_p(s);
  bb = jnrm2_p(z);
  s /=  aa;
  z /=  bb;
  cc = sqrt((1 + dot(s, z)) / 2.0);
  vs = dot(W["v"], s);
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
  }
  W["lambda"] = sqrt(aa * bb) * W["lambda"];
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
std::map<std::string,mat> ntsu_s(std::map<std::string,mat> W, mat s, mat z, int m){
  mat zts, U, V, DiagL, lu;
  vec l;

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
mat sslb_nl(mat s, mat lambda, bool invers){
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
mat sslb_p(mat s, mat lambda, bool invers){
  int n = s.n_rows;
  double a, cc, lx, s0;

  a = jnrm2_p(lambda);
  if(invers){
    lx = dot(lambda, s) / a;
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
mat sslb_s(mat s, mat lambda, bool invers, int m){
  vec ld, ls;

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
mat ssnt_n(mat s, std::map<std::string,mat> W, bool invers){
  mat w;
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
mat ssnt_l(mat s, std::map<std::string,mat> W, bool invers){
  mat w;
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
mat ssnt_p(mat s, std::map<std::string,mat> W, bool invers){
  double a;
  mat w;

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
mat ssnt_s(mat s, std::map<std::string,mat> W, bool invers, bool transp){
  int n, m;
  bool tt;
  mat w, S, a, ans;

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
/*
 * Evaluating a R functions
*/
double feval(mat x, Rcpp::Function Rf){
  double ans;
  ans = Rcpp::as<double>(Rf(Rcpp::wrap(x)));

  return ans;
}
vec geval(mat x, Rcpp::Function Rf){
  vec ans;
  ans = Rcpp::as<vec>(Rf(Rcpp::wrap(x)));

  return ans;
}
mat heval(mat x, Rcpp::Function Rf){
  mat ans;
  ans = Rcpp::as<mat>(Rf(Rcpp::wrap(x)));

  return ans;
}
/*
 * Objective, Gradient and Hessian functions for Risk Parity
*/
double rpp_f0(mat x, mat P, mat mrc){
  double ans = (0.5 * x.t() * P * x - dot(mrc, log(x))).at(0, 0);
  return ans;
}
mat rpp_g0(mat x, mat P, mat mrc){
  mat ans = P * x - mrc / x;
  return ans;
}
mat rpp_h0(mat x, mat P, mat mrc){
  mat ans = zeros(P.n_cols, P.n_cols);
  ans = diagmat(mrc / (x % x)); 
  ans += P;  
  return ans;
}

/*
 * Function value, Gradient and Hessian for geometric programs
*/

std::vector<mat> fgp(mat x, mat F, mat g){
  std::vector<mat> ans;
  double ysum, ymax;
  mat y, fval(1, 1), gval(F.n_cols, 1), hval(F.n_cols, F.n_cols), Fisc;

  y = F * x + g;
  ymax = y.max();
  y = exp(y - ymax);
  ysum = norm(y, 1);
  fval(0, 0) = ymax + log(ysum);

  y /= ysum;
  gval = F.t() * y;

  Fisc = sqrt(diagmat(y)) * (F - ones(F.n_rows, 1) * gval.t());
  hval = Fisc.t() * Fisc;

  ans.push_back(fval);
  ans.push_back(gval);
  ans.push_back(hval);

  return ans;
}
