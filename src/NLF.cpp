#include "NLF.h"

// udot-method
double NLFV::udot(const NLFV& z) const{ // for Rcpp Module
  return arma::dot(u, z.u);
}
// uone-method
NLFV* NLFV::uone() const{
  arma::mat a(dims, 1);
  return new NLFV(a.ones(), dims);
}
// uprd-method
NLFV* NLFV::uprd(const NLFV& z) const{ // for Rcpp Module
    return new NLFV(u % z.u, dims);
}
// uinv-method
NLFV* NLFV::uinv(const NLFV& z) const{ // for Rcpp Module
  return new NLFV(u / z.u, dims);
}
// umss-method
double NLFV::umss() const{
  return u.min();
}
// umsa-method
NLFV* NLFV::umsa(double alpha, bool init) const{ // for Rcpp Module
  arma::mat ua = u;
  if(init){
    for(int i = 0; i < dims; i++){
      ua(i, 0) = u(i, 0) + (1 + alpha);
    }
  } else {
    for(int i = 0; i < dims; i++){
      ua(i, 0) = 1.0 + alpha * u(i, 0);
    }
  }
  return new NLFV(ua, dims);
}
// ntsc-method
NLFS* NLFV::ntsc(const NLFV& z) const{ // for Rcpp Module
  arma::mat dnl(dims, 1), dnli(dims, 1);
  NLFV lambda(dims);

  for(int i = 0; i < dims; i++){
    dnl(i, 0) = sqrt(u(i, 0) / z.u(i, 0));
    dnli(i, 0) = sqrt(z.u(i, 0) / u(i, 0));
    lambda.u(i, 0) = sqrt(u(i, 0) * z.u(i, 0));
  }
  return new NLFS(dnl, dnli, lambda);
}

/*
 * Module for nonlinear related variables
*/
RCPP_MODULE(NLF){
  Rcpp::class_<NLFV>( "NLFV" )
    .constructor("default constructor")
    .constructor<int>("sets the NLF-variable and its dimension")
    .constructor<arma::mat, int>("sets the NLF-variable and its dimension")

    .property("u", &NLFV::get_u, &NLFV::set_u, "Getter and setter for u")
    .property("dims", &NLFV::get_dims, &NLFV::set_dims, "Getter and setter for dims")

    .const_method("uone", &NLFV::uone)
    .const_method("udot", &NLFV::udot)
    .const_method("uprd", &NLFV::uprd)
    .const_method("uinv", &NLFV::uinv)
    .const_method("umss", &NLFV::umss)
    .const_method("umsa", &NLFV::umsa)
    .const_method("ntsc", &NLFV::ntsc)
    ;

  Rcpp::class_<NLFS>( "NLFS" )
    .constructor<arma::mat, arma::mat, NLFV>("sets the NLF-variable and its dimension")

    .property("dnl", &NLFS::get_dnl, &NLFS::set_dnl, "Getter and setter for dnl")
    .property("dnli", &NLFS::get_dnli, &NLFS::set_dnli, "Getter and setter for dnli")
    .property("lambda", &NLFS::get_lambda, &NLFS::set_lambda, "Getter and setter for lambda")
;
}
