#include "NLF.h"

// udot-method
double NLFV::udot(const NLFV& z) const{ 
  return arma::dot(u, z.u);
}
// uone-method
NLFV* NLFV::uone() const{
  arma::mat a(dims, 1);
  return new NLFV(a.ones(), dims);
}
// uprd-method
NLFV* NLFV::uprd(const NLFV& z) const{ 
    return new NLFV(u % z.u, dims);
}
// uinv-method
NLFV* NLFV::uinv(const NLFV& z) const{ 
  return new NLFV(u / z.u, dims);
}
// umss-method
double NLFV::umss() const{
  return u.min();
}
// umsa-method
NLFV* NLFV::umsa(double alpha, bool init) const{ 
  arma::mat ua = u;

  if(init){
    ua(arma::span::all, 0) += (1 + alpha);
  } else {
    ua(arma::span::all, 0) *= alpha;
    ua(arma::span::all, 0) += 1.0;
  }

  return new NLFV(ua, dims);
}
// ntsc-method
NLFS* NLFV::ntsc(const NLFV& z) const{ 
  arma::mat dnl(dims, 1), dnli(dims, 1);
  NLFV lambda(dims);

  dnl = sqrt(u / z.u);
  dnli = sqrt(z.u / u);
  lambda.u = sqrt(u % z.u);

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
