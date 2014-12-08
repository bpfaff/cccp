#include "cccp3.h"
/*
 * Methods / functions for cone variables
*/

/*
 * udot: Inner product of variables in S
*/
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
 * uprd: Product of vectors in S
*/
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
