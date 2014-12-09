#include "PSD.h"

// udot-method
double PSDV::udot(const PSDV& z) const{ 
  double ans = 0.0;
  arma::mat sa = u;
  arma::mat za = z.u;

  sa = arma::reshape(sa, dims, dims);
  za = arma::reshape(za, dims, dims);
  ans = arma::dot(sa.diag(), za.diag());
  for(int i = -1; i > -dims; i--){
    ans += 2 * arma::dot(sa.diag(i), za.diag(i));
  } 

  return ans;
}
// uprd-method
PSDV* PSDV::uprd(const PSDV& z) const{ 
  arma::mat ans(dims, dims);
  arma::mat sa = u;
  arma::mat za = z.u;

  sa.reshape(dims, dims);
  za.reshape(dims, dims);

  ans = 0.5 * (sa * za + za * sa);
  ans.reshape(dims * dims, 1);

  return new PSDV(ans, dims);
}
// uinv-method
PSDV* PSDV::uinv(const PSDV& z) const{ 
  arma::mat sa = u;
  arma::mat za = z.u;

  sa.reshape(dims, dims);
  za.reshape(dims, dims);
  arma::mat ans = sa;

  for(int i = 0; i < dims; i++){
    for(int j = 0; j < dims; j++){
      ans.at(i, j) = sa.at(i, j) * 2.0 / (za.at(i, i) + za.at(j, j));
    }
  }
  ans.reshape(dims * dims, 1);

  return new PSDV(ans, dims);
}
// uone-method
PSDV* PSDV::uone() const{ 
  arma::mat ans = arma::eye(dims, dims);
  ans.reshape(dims * dims, 1);

  return new PSDV(ans, dims);
}
// umss-method
double PSDV::umss() const{ 
  arma::vec eval;
  arma::mat evec;
  arma::mat sa = u;

  sa.reshape(dims, dims);
  arma::eig_sym(eval, evec, sa);

  return -eval[0];
}
// umsa-method
PSDV* PSDV::umsa1(double alpha) const{
  arma::mat sa = u;
  sa.reshape(dims, dims);
  sa.diag() = sa.diag() + (1 + alpha);
  sa.reshape(dims * dims, 1);

  return new PSDV(sa, dims);
}
PSDV* PSDV::umsa2(double alpha, const PSDV& lambda, arma::vec sigma) const{
  arma::mat sa = u;
  arma::mat lu = lambda.u;
  sa.reshape(dims, dims);
  lu.reshape(dims, dims);

  for(int i = 0; i < dims; i++){
    sigma(i) = 1 + alpha * sigma(i);
    sa.col(i) = sa.col(i) * sqrt(sigma(i) / lu(i, i)); 
  }
  sa.reshape(dims * dims, 1);

  return new PSDV(sa, dims);
}
// ntsc-method
PSDS* PSDV::ntsc(const PSDV& z) const{
  arma::mat r(dims, dims), rti(dims, dims);
  PSDV lambda(dims); 
  arma::mat sa = u;
  arma::mat za = z.u;
  arma::mat sc, zc, szc, U, V;
  arma::vec l;

  sa.reshape(dims, dims);
  za.reshape(dims, dims);
  arma::chol(sc, sa);
  arma::chol(zc, za);
  szc = zc * sc.t();
  svd(U, l, V, szc);
  r = zc.i() * U * diagmat(sqrt(l));
  rti = zc.t() * U * diagmat(1.0 / sqrt(l));
  lambda.u = diagmat(l);
  lambda.u.reshape(dims * dims, 1);

  return new PSDS(r, rti, lambda);
}

/*
 * Module for positive semidefinite related variables
*/
RCPP_MODULE(PSD){
  Rcpp::class_<PSDV>( "PSDV" )
    .constructor("default constructor")
    .constructor<int>("sets the PSD-variable and its dimension")
    .constructor<arma::mat, int>("sets the PSD-variable and its dimension")

    .property("u", &PSDV::get_u, &PSDV::set_u, "Getter and setter for u")
    .property("dims", &PSDV::get_dims, &PSDV::set_dims, "Getter and setter for dims")

    .const_method("udot", &PSDV::udot)
    .const_method("uprd", &PSDV::uprd)
    .const_method("uinv", &PSDV::uinv)
    .const_method("uone", &PSDV::uone)
    .const_method("umss", &PSDV::umss)
    .const_method("umsa1", &PSDV::umsa1)
    .const_method("umsa2", &PSDV::umsa2)
    .const_method("ntsc", &PSDV::ntsc)
    ;

  Rcpp::class_<PSDS>( "PSDS" )
    .constructor<arma::mat, arma::mat, PSDV>("sets the PSD-variable and its dimension")

    .property("r", &PSDS::get_r, &PSDS::set_r, "Getter and setter for r")
    .property("rti", &PSDS::get_rti, &PSDS::set_rti, "Getter and setter for rti")
    .property("lambda", &PSDS::get_lambda, &PSDS::set_lambda, "Getter and setter for lambda")
    ;

}
