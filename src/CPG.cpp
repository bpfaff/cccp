#include "CPG.h"

/*
 * Module for control options of optimization routines
*/
RCPP_MODULE(CPG){
  Rcpp::class_<PDV>( "PDV" )
    .constructor("Default constructor")
    .constructor<arma::vec, arma::vec, Rcpp::List, Rcpp::List, double, double>("sets the PDV-values")

    .property("x", &PDV::get_x, &PDV::set_x, "Getter and setter for x")
    .property("y", &PDV::get_y, &PDV::set_y, "Getter and setter for y")
    .property("s", &PDV::get_s, &PDV::set_s, "Getter and setter for s")
    .property("z", &PDV::get_z, &PDV::set_z, "Getter and setter for z")
    .property("kappa", &PDV::get_kappa, &PDV::set_kappa, "Getter and setter for kappa")
    .property("tau", &PDV::get_tau, &PDV::set_tau, "Getter and setter for tau")
    ;
}
