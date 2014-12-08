#include "cccp3.h"
/*
 * Module for SOC related variables
*/
RCPP_MODULE(SOC){
  Rcpp::class_<SOCV>( "SOCV" )
    .constructor("default constructor")
    .constructor<int>("sets the SOC-variable to one-element")
    .constructor<arma::mat, int>("sets the SOC-variable and its dimension")

    .field("u", &SOCV::u)
    .field("dims", &SOCV::dims)

    .method("uone", &uone_s)
    .method("jdot", &jdot_s)
    .method("udot", &udot_s)
    .method("uprd", &uprd_s)
    .method("uinv", &uinv_s)
    .method("umss", &umss_s)
    .method("umsa", &umsa_s)
    .method("ntsc", &ntsc_s)
;

  Rcpp::class_<SOCS>( "SOCS" )
    .constructor<arma::mat, double, SOCV>("sets the SOC-variable and its dimension")

    .field("v", &SOCS::v)
    .field("beta", &SOCS::beta)
    .field("lambda", &SOCS::lambda)
;
}
/*
 * Module for PSD related variables
*/
RCPP_MODULE(PSD){
  Rcpp::class_<PSDV>( "PSDV" )
    .constructor("default constructor")
    .constructor<int>("sets the PSD-variable to one-element")
    .constructor<arma::mat, int>("sets the PSD-variable and its dimension")

    .field("u", &PSDV::u)
    .field("dims", &PSDV::dims)

    .method("uone", &uone_p)
    .method("udot", &udot_p)
    .method("uprd", &uprd_p)
    .method("uinv", &uinv_p)
    .method("umss", &umss_p)
    .method("umsa1", &umsa_p1)
    .method("umsa2", &umsa_p2)
    ;
}
