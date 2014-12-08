#include "cccp3.h"
/*
 * Methods / functions for cone variables
*/

/*
 * udot: Inner product of variables in S
*/
// SOC-variables
double udot_s(SOCV* s, SOCV* z){ // for Rcpp Module
  return arma::dot(s->u, z->u);
}
// PSD-variables
double udot_p(PSDV* s, PSDV* z){ // for Rcpp Module
    double ans = 0.0;
    arma::mat sa = s->u;
    arma::mat za = z->u;

    sa = arma::reshape(sa, s->dims, s->dims);
    za = arma::reshape(za, z->dims, z->dims);
    ans = arma::dot(sa.diag(), za.diag());
    for(int i = -1; i > -s->dims; i--){
      ans += 2 * arma::dot(sa.diag(i), za.diag(i));
    } 

    return ans;
}
/*
 * jdot: Hyperbolic householder transformation for SOC variables
*/
double jdot_s(SOCV* s, SOCV* z){ // for Rcpp Module
  double a = 0.0;

  a = s->u(0, 0) * z->u(0, 0);
  for(int i = 1; i < s->dims; i++){
    a -= s->u(i, 0) * z->u(i, 0);
  }
  return a;
}
/*
 * uprd: Product of vectors in S
*/
// SOC-variables
SOCV uprd_s(SOCV* s, SOCV* z){ // for Rcpp Module
  arma::mat a(s->dims, 1);

  a(0, 0) = arma::dot(s->u, z->u);
  for(int i = 1; i < s->dims; i++){
    a(i, 0) = s->u(0, 0) * z->u(i, 0) + z->u(0, 0) * s->u(i, 0);
  }
  a.reshape(s->dims, 1);

  return SOCV(a, s->dims);
}
// PSD-variables
PSDV uprd_p(PSDV* s, PSDV* z){ // for Rcpp Module
  int dims = s->dims;
  arma::mat a(dims, dims);
  s->u.reshape(dims, dims);
  z->u.reshape(dims, dims);

  a = 0.5 * (s->u * z->u + z->u * s->u);
  a.reshape(dims * dims, 1);

  return PSDV(a, dims);
}
/*
 * uinv: Inverse product of vectors in S
*/
// SOC-variables
SOCV uinv_s(SOCV* s, SOCV* z){ // for Rcpp Module
  int m = s->dims;
  arma::mat a(m, 1);

  double aa = jdot_s(z, z);
  double cc = s->u(0, 0);
  double dd = 0.0;

  for(int i = 1; i < m; i++){
    dd += s->u(i, 0) * z->u(i, 0);
  }
  a(0, 0) = cc * z->u(0, 0) - dd;
  for(int i = 1; i < m; i++){
    a(i, 0) = aa / z->u(0, 0) * s->u(i, 0);
    a(i, 0) += (dd / z->u(0, 0) - cc) * z->u(i, 0);
  }
  for(int i = 0; i < m; i++){
    a(i, 0) /= aa;
  }

  return SOCV(a, m);
}
// PSD-variables
PSDV uinv_p(PSDV* s, PSDV* z){ // for Rcpp Module
  int m = s->dims;
  s->u.reshape(m, m);
  z->u.reshape(z->dims, z->dims);
  arma::mat a = s->u;

  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      a.at(i, j) = s->u.at(i, j) * 2.0 / (z->u.at(i, i) + z->u.at(j, j));
    }
  }
  a.reshape(m * m, 1);

  return PSDV(a, m);
}
/*
 * uone: One-element of vectors in S
*/
// SOC-variables
SOCV uone_s(SOCV* s){ // for Rcpp Module
  arma::mat a(s->dims, 1);
  a = a.zeros();
  a(0, 0) = 1.0;
  return SOCV(a, s->dims);
}
// PSD-variables
PSDV uone_p(PSDV* s){ // for Rcpp Module
  int m = s->dims;
  arma::mat a = arma::eye(m, m);
  a.reshape(m * m, 1);

  return PSDV(a, m);
}
/*
 * umss: Computing maximum step-length for vectors in S
*/
// SOC-variables
SEXP umss_s(SOCV* s){
  double ms = 0.0;

  for(int i = 1; i < s->dims; i++){
    ms += s->u(i, 0) * s->u(i, 0);
  }
  ms = sqrt(ms);
  ms = ms - s->u(0, 0);

  return(Rcpp::List::create(Rcpp::Named("ms") = ms,
			    Rcpp::Named("evd") = R_NilValue));
}
// PSD-variables
SEXP umss_p(PSDV* s){

  arma::vec eval;
  arma::mat evec;
  s->u.reshape(s->dims, s->dims);
  arma::eig_sym(eval, evec, s->u);
  Rcpp::List evd;
  evd = Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(eval),
			   Rcpp::Named("vectors") = Rcpp::wrap(evec));


  return(Rcpp::List::create(Rcpp::Named("ms") = Rcpp::wrap(-eval[0]),
			    Rcpp::Named("evd") = evd));
}
/*
 * umsa: Applying maximum step-size to vectors in S
*/
// SOC-variables
SOCV umsa_s(SOCV* s, double alpha, bool init){
  if(init){
    s->u(0, 0) += 1.0 + alpha;
  } else {
    for(int i = 0; i < s->dims; i++){
      s->u(i, 0) = alpha * s->u(i, 0);
    }
    s->u(0, 0) += 1.0;
  }
  return *s;
}
// PSD-variables (for initial max step length)
PSDV umsa_p1(PSDV* s, double alpha){ 
  s->u.reshape(s->dims, s->dims);
  s->u.diag() = s->u.diag() + (1 + alpha);
  s->u.reshape(s->dims * s->dims, 1);

  return *s;
}
// PSD-variables (for max step length)
PSDV umsa_p2(PSDV* s, double alpha, arma::vec sigma, PSDV* lmbda){
  s->u.reshape(s->dims, s->dims);
  lmbda->u.reshape(lmbda->dims, lmbda->dims);
  for(int i = 0; i < s->dims; i++){
    sigma(i) = 1 + alpha * sigma(i);
    s->u.col(i) = s->u.col(i) * sqrt(sigma(i) / lmbda->u(i, i)); 
  }
  s->u.reshape(s->dims * s->dims, 1);
  
  return *s;
}
/*
 * ntsc: Compute initial NT-scalings and Lagrange Multipliers
*/
// SOC-variables
SOCS ntsc_s(SOCV* s, SOCV* z){
  int m = s->dims;
  arma::mat v(m, 1);
  SOCV lambda(m);
  double aa, bb, cc, dd, beta, szdot;

  aa = sqrt(jdot_s(s, s));
  bb = sqrt(jdot_s(z, z));
  beta = sqrt(aa / bb);
  szdot = udot_s(s, z);
  cc = sqrt((szdot / aa / bb + 1.0) / 2.0);
  v = -z->u / bb;
  v(0, 0) = -v(0, 0);
  for(int i = 0; i < m; i++){
    v(i, 0) += s->u(i, 0) / aa;
    v(i, 0) *= (1.0 / 2.0 / cc);
  }
  v(0, 0) = v(0, 0) + 1.0;
  v *= 1.0 / sqrt(2.0 * v(0, 0));
  lambda.u(0, 0) = cc;
  dd = 2 * cc + s->u(0, 0) / aa + z->u(0, 0) / bb;
  for(int i = 1; i < m; i++){
    lambda.u(i, 0) = s->u(i, 0);
    lambda.u(i, 0) = (cc + z->u(0, 0) / bb) / dd / aa * lambda.u(i, 0);
    lambda.u(i, 0) = (cc + s->u(0, 0) / aa) / dd / bb * z->u(i, 0) + lambda.u(i, 0); 
  }
  lambda.u *= sqrt(aa * bb); 

  return SOCS(v, beta, lambda);
}

