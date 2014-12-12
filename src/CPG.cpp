#include "CPG.h"


// primal objective for QPs
double DQP::pobj(PDV& pdv){
  double ans;
  arma::mat term1;

  term1 = (0.5 * pdv.get_x().t() * P * pdv.get_x());
  ans = term1(0,0) + arma::dot(pdv.get_x(), q);

  return ans;
}
// dual objective for QPs
double DQP::dobj(PDV& pdv){
  double ans;
  arma::mat term1;
  term1(0,0) = 0.0;

  // dobj term for equality constraints
  term1 = pdv.get_y().t() * (A * pdv.get_x() - b);

  // dobj term for inequality constraints
  // needs to be included

  ans = pobj(pdv) + term1(0,0);

  return ans;
}

CPS* DQP::cps(const CTRL& ctrl){
  // Initialising object
  PDV pdv;
  CPS* cps = new CPS();
  Rcpp::NumericVector state = cps->get_state();
  // Case 1: Unconstrained QP
  if((cList.get_K() == 0) && (A.n_rows == 0)){
    pdv.set_x(solve(P, -q));
    cps->set_pdv(pdv);
    state["pobj"] = pobj(pdv);
    cps->set_state(state);
    cps->set_status("optimal");
  }
  return cps;
}

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
    .constructor<arma::vec, arma::vec, std::vector<arma::mat>, std::vector<arma::mat>, double, double>("sets the PDV-values")

    .property("x", &PDV::get_x, &PDV::set_x, "Getter and setter for x")
    .property("y", &PDV::get_y, &PDV::set_y, "Getter and setter for y")
    .property("s", &PDV::get_s, &PDV::set_s, "Getter and setter for s")
    .property("z", &PDV::get_z, &PDV::set_z, "Getter and setter for z")
    .property("kappa", &PDV::get_kappa, &PDV::set_kappa, "Getter and setter for kappa")
    .property("tau", &PDV::get_tau, &PDV::set_tau, "Getter and setter for tau")
    ;


  Rcpp::class_<CONEC>( "CONEC" )
    .constructor("Default constructor")
    .constructor<Rcpp::CharacterVector, Rcpp::List, Rcpp::List, int>("sets the PDV-values")

    .property("conTypes", &CONEC::get_conTypes, &CONEC::set_conTypes, "Getter and setter for conTypes")
    .property("Gmats", &CONEC::get_Gmats, &CONEC::set_Gmats, "Getter and setter for Gmats")
    .property("hvecs", &CONEC::get_hvecs, &CONEC::set_hvecs, "Getter and setter for hvecs")
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
    .method("pobj", &DQP::pobj)
    .method("dobj", &DQP::dobj)
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


