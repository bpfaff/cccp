#include "CPG.h"

/*
 *
 * Methods for Quadratic Programs
 *
*/


/*
Primal objective
*/
double DQP::pobj(PDV& pdv){
  double ans;
  arma::mat term1;

  term1 = (0.5 * pdv.get_x().t() * P * pdv.get_x());
  ans = term1(0,0) + arma::dot(pdv.get_x(), q);

  return ans;
}
/*
Dual objective
*/
double DQP::dobj(PDV& pdv){
  double ans;
  arma::mat term1, term2;
  term1(0,0) = 0.0;
  term2(0,0) = 0.0;

  // dobj term for equality constraints
  term1 = pdv.get_y().t() * (A * pdv.get_x() - b);

  // dobj term for inequality constraints
  if(cList.get_K() > 0){
    for(int i = 0; i < cList.get_K(); i++){
      term2 = term2 + pdv.get_z()[i].t() *			\
	(cList.get_Gmats()[i] * pdv.get_x() - cList.get_hvecs()[i]);
    } 
  }

  ans = pobj(pdv) + term1(0,0) + term2(0,0);

  return ans;
}

/*
  Main routine for solving a Quadratic Program
*/
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
