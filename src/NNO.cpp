#include "NNO.h"

// udot-method
double NNOV::udot(const NNOV& z) const{ 
  return arma::dot(u, z.u);
}
// uone-method
NNOV* NNOV::uone() const{
  arma::mat a(dims, 1);
  return new NNOV(a.ones(), dims);
}
// uprd-method
NNOV* NNOV::uprd(const NNOV& z) const{ 
    return new NNOV(u % z.u, dims);
}
// uinv-method
NNOV* NNOV::uinv(const NNOV& z) const{ 
  return new NNOV(u / z.u, dims);
}
// umss-method
double NNOV::umss() const{
  return u.min();
}
// umsa-method
NNOV* NNOV::umsa(double alpha, bool init) const{ 
  arma::mat ua = u;

  if(init){
    ua(arma::span::all, 0) += (1 + alpha);
  } else {
    ua(arma::span::all, 0) *= alpha;
    ua(arma::span::all, 0) += 1.0;
  }

  return new NNOV(ua, dims);
}
// ntsc-method
NNOS* NNOV::ntsc(const NNOV& z) const{ 
  arma::mat d(dims, 1), di(dims, 1);
  NNOV lambda(dims);

  d = sqrt(u / z.u);
  di = sqrt(z.u / u);
  lambda.u = sqrt(u % z.u);

  return new NNOS(d, di, lambda);
}

/*
 * Module for nonlinear related variables
*/
RCPP_MODULE(NNO){
  Rcpp::class_<NNOV>( "NNOV" )
    .constructor("default constructor")
    .constructor<int>("sets the NNO-variable and its dimension")
    .constructor<arma::mat, int>("sets the NNO-variable and its dimension")

    .property("u", &NNOV::get_u, &NNOV::set_u, "Getter and setter for u")
    .property("dims", &NNOV::get_dims, &NNOV::set_dims, "Getter and setter for dims")

    .const_method("uone", &NNOV::uone)
    .const_method("udot", &NNOV::udot)
    .const_method("uprd", &NNOV::uprd)
    .const_method("uinv", &NNOV::uinv)
    .const_method("umss", &NNOV::umss)
    .const_method("umsa", &NNOV::umsa)
    .const_method("ntsc", &NNOV::ntsc)
    ;

  Rcpp::class_<NNOS>( "NNOS" )
    .constructor<arma::mat, arma::mat, NNOV>("sets the NNO-variable and its dimension")

    .property("d", &NNOS::get_d, &NNOS::set_d, "Getter and setter for d")
    .property("di", &NNOS::get_di, &NNOS::set_di, "Getter and setter for di")
    .property("lambda", &NNOS::get_lambda, &NNOS::set_lambda, "Getter and setter for lambda")
;
}
