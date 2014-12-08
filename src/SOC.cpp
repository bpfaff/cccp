#include "SOC.h"

// udot-method
double SOCV::udot(const SOCV& z) const{ 
  return arma::dot(u, z.u);
}
// jdot-method
double SOCV::jdot(const SOCV& z) const{ 
  double a = 0.0;

  a = u(0, 0) * z.u(0, 0);
  for(int i = 1; i < dims; i++){
    a -= u(i, 0) * z.u(i, 0);
  }

  return a;
}
// uone-method
SOCV* SOCV::uone() const{
  arma::mat a(dims, 1);

  a = a.zeros();
  a(0, 0) = 1.0;

  return new SOCV(a, dims);
}
// uprd-method
SOCV* SOCV::uprd(const SOCV& z) const{ 
  arma::mat a(dims, 1);

  a(0, 0) = arma::dot(u, z.u);
  for(int i = 1; i < dims; i++){
    a(i, 0) = u(0, 0) * z.u(i, 0) + z.u(0, 0) * u(i, 0);
  }
  a.reshape(dims, 1);

  return new SOCV(a, dims);
}
// uinv-method
SOCV* SOCV::uinv(const SOCV& z) const{ 
  arma::mat a(dims, 1);
  double aa = z.jdot(z);
  double cc = u(0, 0);
  double dd = 0.0;

  for(int i = 1; i < dims; i++){
    dd += u(i, 0) * z.u(i, 0);
  }
  a(0, 0) = cc * z.u(0, 0) - dd;
  for(int i = 1; i < dims; i++){
    a(i, 0) = aa / z.u(0, 0) * u(i, 0);
    a(i, 0) += (dd / z.u(0, 0) - cc) * z.u(i, 0);
  }
  for(int i = 0; i < dims; i++){
    a(i, 0) /= aa;
  }

  return new SOCV(a, dims);
}
// umss-method
double SOCV::umss() const{
  double ms = 0.0;

  for(int i = 1; i < dims; i++){
    ms += u(i, 0) * u(i, 0);
  }
  ms = sqrt(ms);
  ms = ms - u(0, 0);

  return ms;
}
// umsa-method
SOCV* SOCV::umsa(double alpha, bool init) const{ 
  arma::mat ua = u;

  if(init){
    ua(arma::span::all, 0) += (1 + alpha);
  } else {
    for(int i = 0; i < dims; i++){
      ua(i, 0) = alpha * ua(i, 0);
    }
    ua(0, 0) += 1.0;
  }

  return new SOCV(ua, dims);
}
// ntsc-method
SOCS* SOCV::ntsc(const SOCV& z) const{ 
  arma::mat v(dims, 1);
  SOCV lambda(dims);
  double aa, bb, cc, dd, beta, szdot;

  aa = sqrt(jdot(*this));
  bb = sqrt(z.jdot(z));
  beta = sqrt(aa / bb);
  szdot = udot(z);
  cc = sqrt((szdot / aa / bb + 1.0) / 2.0);
  v = -z.u / bb;
  v(0, 0) = -v(0, 0);
  for(int i = 0; i < dims; i++){
    v(i, 0) += u(i, 0) / aa;
    v(i, 0) *= (1.0 / 2.0 / cc);
  }
  v(0, 0) = v(0, 0) + 1.0;
  v *= 1.0 / sqrt(2.0 * v(0, 0));
  lambda.u(0, 0) = cc;
  dd = 2 * cc + u(0, 0) / aa + z.u(0, 0) / bb;
  for(int i = 1; i < dims; i++){
    lambda.u(i, 0) = u(i, 0);
    lambda.u(i, 0) = (cc + z.u(0, 0) / bb) / dd / aa * lambda.u(i, 0);
    lambda.u(i, 0) = (cc + u(0, 0) / aa) / dd / bb * z.u(i, 0) + lambda.u(i, 0); 
  }
  lambda.u *= sqrt(aa * bb); 

  return new SOCS(v, beta, lambda);
}
/*
 * Module for second-order cone related variables
*/
RCPP_MODULE(SOC){
  Rcpp::class_<SOCV>( "SOCV" )
    .constructor("default constructor")
    .constructor<int>("sets the SOC-variable and its dimension")
    .constructor<arma::mat, int>("sets the SOC-variable and its dimension")

    .property("u", &SOCV::get_u, &SOCV::set_u, "Getter and setter for u")
    .property("dims", &SOCV::get_dims, &SOCV::set_dims, "Getter and setter for dims")

    .const_method("uone", &SOCV::uone)
    .const_method("udot", &SOCV::udot)
    .const_method("jdot", &SOCV::jdot)
    .const_method("uprd", &SOCV::uprd)
    .const_method("uinv", &SOCV::uinv)
    .const_method("umss", &SOCV::umss)
    .const_method("umsa", &SOCV::umsa)
    .const_method("ntsc", &SOCV::ntsc)
    ;

  Rcpp::class_<SOCS>( "SOCS" )
    .constructor<arma::mat, double, SOCV>("sets the SOC-variable and its dimension")

    .property("v", &SOCS::get_v, &SOCS::set_v, "Getter and setter for v")
    .property("beta", &SOCS::get_beta, &SOCS::set_beta, "Getter and setter for beta")
    .property("lambda", &SOCS::get_lambda, &SOCS::set_lambda, "Getter and setter for lambda")
;
}
