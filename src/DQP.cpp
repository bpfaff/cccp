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
Primal Residuals
*/
arma::mat DQP::rprim(PDV& pdv){
  int p = A.n_rows;
  arma::mat ans(p,1);
  ans.zeros();

  ans = b - A * pdv.get_x();

  return ans;
}

/*
Centrality Resdiuals
*/
std::vector<arma::mat> DQP::rcent(PDV& pdv){
  std::vector<arma::mat> ans;

  for(int i = 0; i < cList.get_K(); i++){
    ans[i] = pdv.get_s()[i] + cList.get_Gmats()[i] * pdv.get_x() - cList.get_hvecs()[i];
  }

  return ans;
}

/*
Dual Residuals
*/
arma::mat DQP::rdual(PDV& pdv){
  int n = P.n_rows;
  arma::mat Gz(n,1);
  arma::mat Ay(n,1);
  arma::mat ans(n,1);
  Gz.zeros();
  Ay.zeros();
  ans.zeros();

  if(cList.get_K() > 0){
    for(int i = 0; i < cList.get_K(); i++){
      Gz = Gz + cList.get_Gmats()[i].t() * pdv.get_z()[i];
    } 
  }
  if(A.n_rows > 0){
    Ay = A.t() * pdv.get_y();
  }
  ans = P * pdv.get_x() + q + Gz + Ay;

  return ans;
}
/*
Certificate of primal infeasibilty
*/
double DQP::certp(PDV& pdv){
  double nomin, denom, ans1 = 0.0, ans2 = 0.0, ans = 0.0;

  nomin = arma::norm(rprim(pdv));
  denom = std::max(1.0, arma::norm(b));
  ans1 = nomin / denom;

  if(cList.get_K() > 0){
    std::vector<arma::mat> rz;
    rz = rcent(pdv);
    nomin = 0.0;
    denom = std::max(1.0, arma::norm(q));
    for(int i = 0; i < cList.get_K(); i++){
      nomin += arma::norm(rz[i]);
    } 
    ans2 = nomin / denom;
  }
  ans = std::max(ans1, ans2);

  return ans;
}

/*
Certificate of dual infeasibilty
*/
double DQP::certd(PDV& pdv){
  double nomin, denom, ans;

  int n = P.n_rows;
  arma::mat Gz(n,1);
  Gz.zeros();

  if(cList.get_K() > 0){
    for(int i = 0; i < cList.get_K(); i++){
      Gz = Gz + cList.get_Gmats()[i].t() * pdv.get_z()[i];
    } 
  }

  nomin = arma::norm(P * pdv.get_x() + Gz + A.t() * pdv.get_y() + q);
  denom = std::max(1.0, arma::norm(q));
  ans = nomin / denom;

  return ans;
}

/*
  Main routine for solving a Quadratic Program
*/
CPS* DQP::cps(CTRL& ctrl){
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
  // Case 2: Equality constrained QP
  if((cList.get_K() == 0) && (A.n_rows > 0)){
    arma::mat Pi, PiA, Piq, S, Si;
    double ftol = ctrl.get_feastol();
    try{
      Pi = inv(P);
    } catch(std::runtime_error &ex){
      forward_exception_to_r(ex);
    } catch(...){
      ::Rf_error("C++ exception (unknown reason)");
    }
    Piq = Pi * q;
    PiA = Pi * A.t();
    S = -A * PiA;
    try{
      Si = inv(S);
    } catch(std::runtime_error &ex){
      throw std::range_error("Inversion of Schur complement failed.");
      forward_exception_to_r(ex);
    } catch(...) {
      ::Rf_error("C++ exception (unknown reason)");
    }
    pdv.set_y(Si * (A * Piq + b));
    pdv.set_x(Pi * (-A.t() * pdv.get_y() - q));
    cps->set_pdv(pdv);
    state["pobj"] = pobj(pdv);
    state["dobj"] = dobj(pdv);
    state["certp"] = certp(pdv);
    state["certd"] = certd(pdv);
    cps->set_state(state);
    if((state["certp"] <= ftol) && (state["certd"] <= ftol)){
      cps->set_status("optimal");
    } else {
      cps->set_status("optimal");
    }
  }
  return cps;
}
