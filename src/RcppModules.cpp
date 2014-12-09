/*
#include "cccp3.h"

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
*/
