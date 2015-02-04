#include "CPG.h"
using namespace arma;
/*
 * Module for control options of optimization routines
*/
RCPP_MODULE(CPG){

  Rcpp::class_<CTRL>( "CTRL" )
    .constructor("Default constructor")
    .constructor<Rcpp::List>("sets the CTRL-options")

    .property("params", &CTRL::get_params, &CTRL::set_params, "Getter and setter for control parameters")
    ;

  Rcpp::class_<PDV>( "PDV" )
    .constructor("Default constructor")
    .constructor<mat, mat, mat, mat, double, double>("PDV-values; stacked slack variables.")

    .property("x", &PDV::get_x, &PDV::set_x, "Getter and setter for x")
    .property("y", &PDV::get_y, &PDV::set_y, "Getter and setter for y")
    .property("s", &PDV::get_s, &PDV::set_s, "Getter and setter for s")
    .property("z", &PDV::get_z, &PDV::set_z, "Getter and setter for z")
    .property("kappa", &PDV::get_kappa, &PDV::set_kappa, "Getter and setter for kappa")
    .property("tau", &PDV::get_tau, &PDV::set_tau, "Getter and setter for tau")
    ;

  Rcpp::class_<CONEC>( "CONEC" )
    .constructor("Default constructor")
    .constructor<std::vector<std::string>, mat, mat, umat, uvec, int, int>("cone constraints")
    .constructor<int>("no cone constraints, setting n")

    .property("cone", &CONEC::get_cone, &CONEC::set_cone, "Getter and setter for cone types")
    .property("G", &CONEC::get_G, &CONEC::set_G, "Getter and setter for Gmats")
    .property("h", &CONEC::get_h, &CONEC::set_h, "Getter and setter for hmats")
    .property("sidx", &CONEC::get_sidx, &CONEC::set_sidx, "Getter and setter for sidx")
    .property("dims", &CONEC::get_dims, &CONEC::set_dims, "Getter and setter for dims")
    .property("K", &CONEC::get_K, &CONEC::set_K, "Getter and setter for K")
    .property("n", &CONEC::get_n, &CONEC::set_n, "Getter and setter for K")
 
    .method("initpdv", &CONEC::initpdv)

   ;

  Rcpp::class_<DCP>( "DCP" )
    .constructor("Default constructor")
    .constructor<mat, CONEC, Rcpp::List, mat, mat>("sets the DCP-values")

    .property("x0", &DCP::get_x0, &DCP::set_x0, "Getter and setter for x0")
    .property("cList", &DCP::get_cList, &DCP::set_cList, "Getter and setter for cList") 
    .property("nList", &DCP::get_nList, &DCP::set_nList, "Getter and setter for nList")
    .property("A", &DCP::get_A, &DCP::set_A, "Getter and setter for A")
    .property("b", &DCP::get_b, &DCP::set_b, "Getter and setter for b")

    .method("cps", &DCP::cps)
    ;

  Rcpp::class_<DLP>( "DLP" )
    .constructor("Default constructor")
    .constructor<vec, mat, vec, CONEC>("sets the DLP-values")

    .property("q", &DLP::get_q, &DLP::set_q, "Getter and setter for q")
    .property("A", &DLP::get_A, &DLP::set_A, "Getter and setter for A")
    .property("b", &DLP::get_b, &DLP::set_b, "Getter and setter for b")
    .property("cList", &DLP::get_cList, &DLP::set_cList, "Getter and setter for cList")

    .method("cps", &DLP::cps)
    ;

  Rcpp::class_<DNL>( "DNL" )
    .constructor("Default constructor")
    .constructor<vec, mat, vec, CONEC, mat, Rcpp::List>("sets the DNL-values")

    .property("q", &DNL::get_q, &DNL::set_q, "Getter and setter for q")
    .property("A", &DNL::get_A, &DNL::set_A, "Getter and setter for A")
    .property("b", &DNL::get_b, &DNL::set_b, "Getter and setter for b")
    .property("cList", &DNL::get_cList, &DNL::set_cList, "Getter and setter for cList")
    .property("x0", &DNL::get_x0, &DNL::set_x0, "Getter and setter for x0")
    .property("nList", &DNL::get_nList, &DNL::set_nList, "Getter and setter for nList")

    .method("cps", &DNL::cps)
    ;

  Rcpp::class_<DQP>( "DQP" )
    .constructor("Default constructor")
    .constructor<mat, vec, mat, vec, CONEC>("sets the DQP-values")

    .property("P", &DQP::get_P, &DQP::set_P, "Getter and setter for P")
    .property("q", &DQP::get_q, &DQP::set_q, "Getter and setter for q")
    .property("A", &DQP::get_A, &DQP::set_A, "Getter and setter for A")
    .property("b", &DQP::get_b, &DQP::set_b, "Getter and setter for b")
    .property("cList", &DQP::get_cList, &DQP::set_cList, "Getter and setter for cList")

    .method("cps", &DQP::cps)
    ;

  Rcpp::class_<CPS>( "CPS" )
    .constructor("Default constructor")
    .constructor<PDV, Rcpp::NumericVector, Rcpp::String, int, umat>("sets the CPS-values")

    .property("pdv", &CPS::get_pdv, &CPS::set_pdv, "Getter and setter for pdv")
    .property("state", &CPS::get_state, &CPS::set_state, "Getter and setter for state")
    .property("status", &CPS::get_status, &CPS::set_status, "Getter and setter for status")
    .property("niter", &CPS::get_niter, &CPS::set_niter, "Getter and setter for niter")
    .property("sidx", &CPS::get_sidx, &CPS::set_sidx, "Getter and setter for sidx")
    ;

  Rcpp::function("rpp", &rpp);
  Rcpp::function("gpp", &gpp);
}


