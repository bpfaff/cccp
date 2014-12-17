#include "CPG.h"

/*
 * Module for control options of optimization routines
*/
RCPP_MODULE(CPG){
  Rcpp::class_<CTRL>( "CTRL" )
    .constructor<int, double, double, double, bool>("sets the CTRL-options")

    .property("maxiters", &CTRL::get_maxiters, &CTRL::set_maxiters, "Getter and setter for maxiters")
    .property("abstol", &CTRL::get_abstol, &CTRL::set_abstol, "Getter and setter for abstol")
    .property("reltol", &CTRL::get_reltol, &CTRL::set_reltol, "Getter and setter for reltol")
    .property("feastol", &CTRL::get_feastol, &CTRL::set_feastol, "Getter and setter for feastol")
    .property("trace", &CTRL::get_trace, &CTRL::set_trace, "Getter and setter for trace")
    ;

  Rcpp::class_<PDV>( "PDV" )
    .constructor("Default constructor")
    .constructor<arma::mat, arma::mat, std::vector<arma::mat>, std::vector<arma::mat>, double, double>("sets the PDV-values")

    .property("x", &PDV::get_x, &PDV::set_x, "Getter and setter for x")
    .property("y", &PDV::get_y, &PDV::set_y, "Getter and setter for y")
    .property("s", &PDV::get_s, &PDV::set_s, "Getter and setter for s")
    .property("z", &PDV::get_z, &PDV::set_z, "Getter and setter for z")
    .property("kappa", &PDV::get_kappa, &PDV::set_kappa, "Getter and setter for kappa")
    .property("tau", &PDV::get_tau, &PDV::set_tau, "Getter and setter for tau")
    ;


  Rcpp::class_<CONEC>( "CONEC" )
    .constructor("Default constructor")
    .constructor<Rcpp::CharacterVector, Rcpp::List, Rcpp::List, Rcpp::IntegerVector, int>("sets the inequality constraints")

    .property("conTypes", &CONEC::get_conTypes, &CONEC::set_conTypes, "Getter and setter for conTypes")
    .property("Gmats", &CONEC::get_Gmats, &CONEC::set_Gmats, "Getter and setter for Gmats")
    .property("hvecs", &CONEC::get_hvecs, &CONEC::set_hvecs, "Getter and setter for hvecs")
    .property("dims", &CONEC::get_dims, &CONEC::set_dims, "Getter and setter for dims")
    .property("K", &CONEC::get_K, &CONEC::set_K, "Getter and setter for hvecs")
    ;

  Rcpp::class_<DQP>( "DQP" )
    .constructor("Default constructor")
    .constructor<arma::mat, arma::vec, arma::mat, arma::vec, CONEC>("sets the DQP-values")

    .property("P", &DQP::get_P, &DQP::set_P, "Getter and setter for P")
    .property("q", &DQP::get_q, &DQP::set_q, "Getter and setter for q")
    .property("A", &DQP::get_A, &DQP::set_A, "Getter and setter for A")
    .property("b", &DQP::get_b, &DQP::set_b, "Getter and setter for b")
    .property("cList", &DQP::get_cList, &DQP::set_cList, "Getter and setter for cList")

    .method("cps", &DQP::cps)
    ;

  Rcpp::class_<CPS>( "CPS" )
    .constructor("Default constructor")
    .constructor<PDV, Rcpp::NumericVector, Rcpp::String, int>("sets the CPS-values")

    .property("pdv", &CPS::get_pdv, &CPS::set_pdv, "Getter and setter for pdv")
    .property("state", &CPS::get_state, &CPS::set_state, "Getter and setter for state")
    .property("status", &CPS::get_status, &CPS::set_status, "Getter and setter for status")
    .property("niter", &CPS::get_niter, &CPS::set_niter, "Getter and setter for niter")
    ;
}


