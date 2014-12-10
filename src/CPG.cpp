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

  Rcpp::class_<DQP>( "DQP" )
    .constructor("Default constructor")
    .constructor<arma::mat, arma::vec, arma::mat, arma::vec, Rcpp::List>("sets the DQP-values")

    .property("P", &DQP::get_P, &DQP::set_P, "Getter and setter for P")
    .property("q", &DQP::get_q, &DQP::set_q, "Getter and setter for q")
    .property("A", &DQP::get_A, &DQP::set_A, "Getter and setter for A")
    .property("b", &DQP::get_b, &DQP::set_b, "Getter and setter for b")
    .property("cList", &DQP::get_cList, &DQP::set_cList, "Getter and setter for cList")
    ;
}
