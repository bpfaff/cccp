#include "cccp3.h"

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
  arma::mat term1(1,1);
  term1(0,0) = 0.0;

  term1 = (0.5 * pdv.x.t() * P * pdv.x);
  ans = term1(0,0) + arma::dot(pdv.x, q);

  return ans;
}
/*
Dual objective
*/
double DQP::dobj(PDV& pdv){
  double ans;
  arma::mat term1(1,1), term2(1,1);
  term1(0,0) = 0.0;
  term2(0,0) = 0.0;

  // dobj term for equality constraints
  term1 = pdv.y.t() * (A * pdv.x - b);
  // dobj term for inequality constraints
  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      term2 = term2 + pdv.z[i].t() *			\
	(cList.Gmats[i] * pdv.x - cList.hmats[i]);
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

  ans = b - A * pdv.x;

  return ans;
}
/*
Centrality Resdiuals
*/
std::vector<arma::mat> DQP::rcent(PDV& pdv){
  std::vector<arma::mat> ans;

  for(int i = 0; i < cList.K; i++){
    ans[i] = pdv.s[i] + cList.Gmats[i] * pdv.x - cList.hmats[i];
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

  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      Gz = Gz + cList.Gmats[i].t() * pdv.z[i];
    } 
  }
  if(A.n_rows > 0){
    Ay = A.t() * pdv.y;
  }
  ans = P * pdv.x + q + Gz + Ay;

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

  if(cList.K > 0){
    std::vector<arma::mat> rz;
    rz = rcent(pdv);
    nomin = 0.0;
    denom = std::max(1.0, arma::norm(q));
    for(int i = 0; i < cList.K; i++){
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

  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      Gz = Gz + cList.Gmats[i].t() * pdv.z[i];
    } 
  }

  nomin = arma::norm(P * pdv.x + Gz + A.t() * pdv.y + q);
  denom = std::max(1.0, arma::norm(q));
  ans = nomin / denom;

  return ans;
}
/*
Initializing PDV
*/
PDV* DQP::initpdv(){
  PDV* pdv = new PDV();
  std::vector<arma::mat> s, z;
  arma::mat ans;
  int n = P.n_cols;

  pdv->x = arma::zeros(n,1);
  pdv->y = arma::zeros(A.n_rows,1);
  for(int i = 0; i < cList.K; i++){
    if((cList.conTypes[i] == "NLFC") || (cList.conTypes[i] == "NNOC")){
      ans = arma::ones(cList.dims[i], 1);
    } else if(cList.conTypes[i] == "SOCC"){
      ans = arma::zeros(cList.dims[i], 1);
      ans.at(0,0) = 1.0;
    } else if(cList.conTypes[i] == "PSDC") {
      ans = arma::eye(cList.dims[i],cList.dims[i]);
      ans.reshape(cList.dims[i] * cList.dims[i], 1);
    } else {
      ans = arma::zeros(cList.dims[i], 1);
    }
    s[i] = ans;
    z[i] = ans;
  }
  pdv->s = s;
  pdv->z = s;
  pdv->tau = 1.0;
  pdv->kappa = 1.0;

  return pdv;
}
/*
Initial Nesterov-Todd scalings
*/
std::vector<std::map<std::string,arma::mat> > DQP::initnts(){
  std::vector<std::map<std::string,arma::mat> > WList;
  std::map<std::string,arma::mat> W;
  arma::mat ans;

  for(int i = 0; i < cList.K; i++){
    if(cList.conTypes[i] == "NLFC"){
      ans = arma::ones(cList.dims[i],1);
      W["dnl"] = ans;
      W["dnli"] = ans;
      ans = arma::zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.conTypes[i] == "NNOC"){
      ans = arma::ones(cList.dims[i],1);
      W["d"] = ans;
      W["di"] = ans;
      ans = arma::zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.conTypes[i] == "SOCC"){
      ans = arma::ones(1,1);
      W["beta"] = ans;
      ans = arma::zeros(cList.dims[i],1);
      ans.at(0,0) = 1.0;
      W["v"] = ans;
      ans = arma::zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.conTypes[i] == "PSDC"){
      ans = arma::eye(cList.dims[i],cList.dims[i]);
      ans.reshape(cList.dims[i] * cList.dims[i], 1);
      W["r"] = ans;
      W["rti"] = ans;
      ans = arma::zeros(cList.dims[i] * cList.dims[i], 1);
      W["lambda"] = ans;
    }
    WList[i] = W;
  }

  return WList;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
arma::mat DQP::gwwg(std::vector<std::map<std::string,arma::mat> > WList){
  int n = P.n_cols;
  arma::mat gwwg(n,n), temp(n,n), witg, wiwitg;
  gwwg.zeros();
  temp.zeros();

  for(int i = 0; i < cList.K; i++){
    if(cList.conTypes[i] == "NLFC"){
      witg = ssnt_n(cList.Gmats[i], WList[i], true);
      wiwitg = ssnt_n(witg, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitg;
    } else if(cList.conTypes[i] == "NNOC"){
      witg = ssnt_l(cList.Gmats[i], WList[i], true);
      wiwitg = ssnt_l(witg, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitg;
    } else if(cList.conTypes[i] == "SOCC"){
      witg = ssnt_p(cList.Gmats[i], WList[i], true);
      wiwitg = ssnt_p(witg, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitg;
    } else if(cList.conTypes[i] == "PSDC"){
      witg = ssnt_s(cList.Gmats[i], WList[i], true, true);
      wiwitg = ssnt_s(witg, WList[i], true, false);
      temp = cList.Gmats[i].t() * wiwitg;
    }
    gwwg = gwwg + temp;
  }

  return gwwg;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
arma::mat DQP::gwwz(std::vector<std::map<std::string,arma::mat> > WList, std::vector<arma::mat> z){
  int n = P.n_cols;
  arma::mat gwwz(n,1), temp(n,1), witz, wiwitz;
  gwwz.zeros();
  temp.zeros();

  for(int i = 0; i < cList.K; i++){
    if(cList.conTypes[i] == "NLFC"){
      witz = ssnt_n(z[i], WList[i], true);
      wiwitz = ssnt_n(witz, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitz;
    } else if(cList.conTypes[i] == "NNOC"){
      witz = ssnt_l(z[i], WList[i], true);
      wiwitz = ssnt_l(witz, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitz;
    } else if(cList.conTypes[i] == "SOCC"){
      witz = ssnt_p(z[i], WList[i], true);
      wiwitz = ssnt_p(witz, WList[i], true);
      temp = cList.Gmats[i].t() * wiwitz;
    } else if(cList.conTypes[i] == "PSDC"){
      witz = ssnt_s(z[i], WList[i], true, true);
      wiwitz = ssnt_s(witz, WList[i], true, false);
      temp = cList.Gmats[i].t() * wiwitz;
    }
    gwwz = gwwz + temp;
  }

  return gwwz;
}
/*
Solving 'KKT-System'
*/

/*
  Main routine for solving a Quadratic Program
*/
CPS* DQP::cps(CTRL& ctrl){
  // Initialising object
  PDV* pdv = initpdv();
  CPS* cps = new CPS();
  Rcpp::NumericVector state = cps->get_state();
  // Case 1: Unconstrained QP
  if((cList.K == 0) && (A.n_rows == 0)){
    pdv->x = solve(P, -q);
    cps->set_pdv(*pdv);
    state["pobj"] = pobj(*pdv);
    cps->set_state(state);
    cps->set_status("optimal");
  }
  // Case 2: Equality constrained QP
  if((cList.K == 0) && (A.n_rows > 0)){
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

    pdv->y = Si * (A * Piq + b);
    pdv->x = Pi * (-(A.t() * pdv->y) - q);
    cps->set_pdv(*pdv);
    state["pobj"] = pobj(*pdv);
    state["dobj"] = dobj(*pdv);
    state["certp"] = certp(*pdv);
    state["certd"] = certd(*pdv);
    cps->set_state(state);
    if((state["certp"] <= ftol) && (state["certd"] <= ftol)){
      cps->set_status("optimal");
    } else {
      cps->set_status("unknown");
    }
  }
  // Case 3: At least inequality constrained QP
  // Initialising state variables
  int m = sum(cList.dims);
  std::map<std::string,double> cvgdvals;
  cvgdvals["pobj"] = NA_REAL;
  cvgdvals["dobj"] = NA_REAL;
  cvgdvals["pinf"] = NA_REAL;
  cvgdvals["dinf"] = NA_REAL;
  cvgdvals["dgap"] = NA_REAL;
  // Initialising LHS and RHS matrices
  int n = P.n_cols;
  int sizeLHS = A.n_rows + A.n_cols;
  arma::mat LHS(sizeLHS, sizeLHS);
  LHS.zeros();
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(n, 0, sizeLHS-1, n-1) = A;
    LHS.submat(0, n, n-1, sizeLHS-1) = A.t();
  }
  arma::mat RHS(sizeLHS, 1);
  // Initialising Nesterov-Todd scalings
  std::vector<std::map<std::string,arma::mat> > WList;
  WList = initnts();
  // Initialising PDV for determining (first) interior point
  pdv->x = -q;
  pdv->y = b;
  pdv->z = cList.hmats;

  Rcpp::Rcout << "test" << std::endl;
 
  return cps;
}
