#include "cccp3.h"
/*
 * Module for nonlinear related variables
*/
RCPP_MODULE(NLF){
  Rcpp::class_<NLFV>( "NLFV" )
    .constructor("default constructor")
    .constructor<int>("sets the NLF-variable to one-element")
    .constructor<arma::mat, int>("sets the NLF-variable and its dimension")

    .field("u", &NLFV::u)
    .field("dims", &NLFV::dims)

    .method("uone", &uone_n)
    .method("udot", &udot_n)
    .method("uprd", &uprd_n)
    .method("uinv", &uinv_n)
    .method("umss", &umss_n)
    .method("umsa", &umsa_n)
    .method("ntsc", &ntsc_n)
;

  Rcpp::class_<NLFS>( "NLFS" )
    .constructor<arma::mat, arma::mat, NLFV>("sets the NLF-variable and its dimension")

    .field("dnl", &NLFS::dnl)
    .field("dnli", &NLFS::dnli)
    .field("lambda", &NLFS::lambda)
;
}
/*
 * Module for NNO related variables
*/
RCPP_MODULE(NNO){
  Rcpp::class_<NNOV>( "NNOV" )
    .constructor("default constructor")
    .constructor<int>("sets the NNO-variable to one-element")
    .constructor<arma::mat, int>("sets the NNO-variable and its dimension")

    .field("u", &NNOV::u)
    .field("dims", &NNOV::dims)

    .method("uone", &uone_l)
    .method("udot", &udot_l)
    .method("uprd", &uprd_l)
    .method("uinv", &uinv_l)
    .method("umss", &umss_l)
    .method("umsa", &umsa_l)
    .method("ntsc", &ntsc_l)
    ;

  Rcpp::class_<NNOS>( "NNOS" )
    .constructor<arma::mat, arma::mat, NNOV>("sets the NNO-variable and its dimension")

    .field("d", &NNOS::d)
    .field("di", &NNOS::di)
    .field("lambda", &NNOS::lambda)
;
}
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
