#include "cccp3.h"

/*
 *
 * Methods for Quadratic Programs
 *
*/
using namespace arma;

/*
Primal objective
*/
double DQP::pobj(PDV& pdv){
  double ans;
  mat term1(1,1);
  term1(0,0) = 0.0;

  term1 = (0.5 * pdv.x.t() * P * pdv.x);
  ans = term1(0,0) + dot(pdv.x, q);

  return ans;
}
/*
Dual objective
*/
double DQP::dobj(PDV& pdv){
  double ans;
  mat term1(1,1), term2(1,1);
  term1(0,0) = 0.0;
  term2(0,0) = 0.0;

  // dobj term for equality constraints
  term1 = pdv.y.t() * (A * pdv.x - b);
  // dobj term for inequality constraints
  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      term2 = term2 + pdv.z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * \
	(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * pdv.x - \
	 cList.h(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all));
    } 
  }
  ans = pobj(pdv) + term1(0,0) + term2(0,0);

  return ans;
}
/*
Primal Residuals
*/
mat DQP::rprim(PDV& pdv){
  int p = A.n_rows;
  mat ans(p,1);
  ans.zeros();

  ans = b - A * pdv.x;

  return ans;
}
/*
Centrality Resdiuals
*/
mat DQP::rcent(PDV& pdv){
  mat ans(cList.G.n_rows, P.n_cols);

  for(int i = 0; i < cList.K; i++){
    ans(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = pdv.s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) + cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * pdv.x - cList.h(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all);
  }

  return ans;
}
/*
Dual Residuals
*/
mat DQP::rdual(PDV& pdv){
  int n = P.n_rows;
  mat Gz(n,1);
  mat Ay(n,1);
  mat ans(n,1);
  Gz.zeros();
  Ay.zeros();
  ans.zeros();

  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      Gz = Gz + cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * pdv.z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all);
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

  nomin = norm(rprim(pdv));
  denom = std::max(1.0, norm(b));
  ans1 = nomin / denom;

  if(cList.K > 0){
    mat rz;
    rz = rcent(pdv);
    nomin = 0.0;
    denom = std::max(1.0, norm(q));
    for(int i = 0; i < cList.K; i++){
      nomin += norm(rz(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all));
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
  mat Gz(n,1);
  Gz.zeros();

  if(cList.K > 0){
    for(int i = 0; i < cList.K; i++){
      Gz = Gz + cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * pdv.z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all);
    } 
  }

  nomin = norm(P * pdv.x + Gz + A.t() * pdv.y + q);
  denom = std::max(1.0, norm(q));
  ans = nomin / denom;

  return ans;
}
/*
Initializing PDV
*/
PDV* DQP::initpdv(){
  PDV* pdv = new PDV();
  mat s(cList.G.n_rows, 1);
  mat ans;
  int n = P.n_cols;

  pdv->x = zeros(n,1);
  pdv->y = zeros(A.n_rows,1);
  for(int i = 0; i < cList.K; i++){
    if((cList.cone[i] == "NLFC") || (cList.cone[i] == "NNOC")){
      ans = ones(cList.dims[i], 1);
    } else if(cList.cone[i] == "SOCC"){
      ans = zeros(cList.dims[i], 1);
      ans.at(0,0) = 1.0;
    } else if(cList.cone[i] == "PSDC") {
      ans = eye(cList.dims[i],cList.dims[i]);
      ans.reshape(cList.dims[i] * cList.dims[i], 1);
    } else {
      ans = zeros(cList.dims[i], 1);
    }
    s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = ans;
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
std::vector<std::map<std::string,mat> > DQP::initnts(){
  std::vector<std::map<std::string,mat> > WList;
  std::map<std::string,mat> W;
  mat ans;

  for(int i = 0; i < cList.K; i++){
    if(cList.cone[i] == "NLFC"){
      ans = ones(cList.dims[i],1);
      W["dnl"] = ans;
      W["dnli"] = ans;
      ans = zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.cone[i] == "NNOC"){
      ans = ones(cList.dims[i],1);
      W["d"] = ans;
      W["di"] = ans;
      ans = zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.cone[i] == "SOCC"){
      ans = ones(1,1);
      W["beta"] = ans;
      ans = zeros(cList.dims[i],1);
      ans.at(0,0) = 1.0;
      W["v"] = ans;
      ans = zeros(cList.dims[i],1);
      W["lambda"] = ans;
    } else if(cList.cone[i] == "PSDC"){
      ans = eye(cList.dims[i],cList.dims[i]);
      ans.reshape(cList.dims[i] * cList.dims[i], 1);
      W["r"] = ans;
      W["rti"] = ans;
      ans = zeros(cList.dims[i] * cList.dims[i], 1);
      W["lambda"] = ans;
    }
    WList[i] = W;
  }

  return WList;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
mat DQP::gwwg(std::vector<std::map<std::string,mat> > WList){
  int n = P.n_cols;
  mat gwwg(n,n), temp(n,n), witg, wiwitg;
  gwwg.zeros();
  temp.zeros();

  for(int i = 0; i < cList.K; i++){
    if(cList.cone[i] == "NLFC"){
      witg = ssnt_n(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_n(witg, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cList.cone[i] == "NNOC"){
      witg = ssnt_l(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_l(witg, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cList.cone[i] == "SOCC"){
      witg = ssnt_p(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_p(witg, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cList.cone[i] == "PSDC"){
      witg = ssnt_s(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true, true);
      wiwitg = ssnt_s(witg, WList[i], true, false);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitg;
    }
    gwwg = gwwg + temp;
  }

  return gwwg;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
mat DQP::gwwz(std::vector<std::map<std::string,mat> > WList, mat z){
  int n = P.n_cols;
  mat gwwz(n,1), temp(n,1), witz, wiwitz;
  gwwz.zeros();
  temp.zeros();

  for(int i = 0; i < cList.K; i++){
    if(cList.cone[i] == "NLFC"){
      witz = ssnt_n(z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_n(witz, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cList.cone[i] == "NNOC"){
      witz = ssnt_l(z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_l(witz, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cList.cone[i] == "SOCC"){
      witz = ssnt_p(z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_p(witz, WList[i], true);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cList.cone[i] == "PSDC"){
      witz = ssnt_s(z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], true, true);
      wiwitz = ssnt_s(witz, WList[i], true, false);
      temp = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all).t() * wiwitz;
    }
    gwwz = gwwz + temp;
  }

  return gwwz;
}
/*
Nesterov-Todd scaling for all inequality constraints
*/
mat DQP::ssnt(mat s, std::vector<std::map<std::string,mat> > WList, bool invers, bool transp){

  for(int i = 0; i < cList.K; i++){
    if(cList.cone[i] == "NLFC"){
      s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = ssnt_n(s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cList.cone[i] == "NNOC"){
      s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = ssnt_l(s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cList.cone[i] == "SOCC"){
      s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = ssnt_p(s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cList.cone[i] == "PSDC"){
      s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = ssnt_s(s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all)(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), WList[i], invers, transp);
    }
  }

  return s;
}
/*
Solving 'KKT-System'
*/
PDV* DQP::sxyz(PDV* pdv, mat LHS, mat RHS, std::vector<std::map<std::string,mat> > WList){
  int n = P.n_cols;
  mat lhs1, rhs1, ans;

  lhs1 = gwwg(WList);
  LHS.submat(0, 0, n-1, n-1) = P + lhs1;
  rhs1 = gwwz(WList, pdv->z);
  RHS.submat(0, 0, n - 1, 0) = pdv->x + rhs1;
  RHS.submat(n, 0, RHS.n_rows - 1, 0) = pdv->y;
  ans = solve(LHS, RHS);
  pdv->x = ans.submat(0, 0, n - 1, 0);
  pdv->y = ans.submat(n, 0, RHS.n_rows - 1, 0);
  for(int i = 0; i < cList.K; i++){
    pdv->z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * pdv->x - pdv->z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all);
  }
  pdv->z = ssnt(pdv->z, WList, true, true);

  return pdv;
}
/*
  Main routine for solving a Quadratic Program
*/
CPS* DQP::cps(CTRL& ctrl){
  // Initialising object
  PDV* pdv = initpdv();
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  Rcpp::NumericVector state = cps->get_state();
  // Case 1: Unconstrained QP
  if((cList.K == 0) && (A.n_rows == 0)){
    pdv->x = solve(P, -q);
    cps->set_pdv(*pdv);
    state["pobj"] = pobj(*pdv);
    cps->set_state(state);
    cps->set_status("optimal");
    return cps;
  }
  // Case 2: Equality constrained QP
  if((cList.K == 0) && (A.n_rows > 0)){
    mat Pi, PiA, Piq, S, Si;
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
    return cps;
  }
  // Case 3: At least inequality constrained QP
    /* 

  if(cList.K > 0){
    // Initialising state variables
    //int m = sum(cList.dims);
   std::map<std::string,double> cvgdvals;
    cvgdvals["pobj"] = NA_REAL;
    cvgdvals["dobj"] = NA_REAL;
    cvgdvals["pinf"] = NA_REAL;
    cvgdvals["dinf"] = NA_REAL;
    cvgdvals["dgap"] = NA_REAL;

    // Initialising LHS and RHS matrices
    int n = P.n_cols;
    int sizeLHS = A.n_rows + A.n_cols;
    mat LHS(sizeLHS, sizeLHS);
    LHS.zeros();
    if(A.n_rows > 0){ // equality constraints
      LHS.submat(n, 0, sizeLHS-1, n-1) = A;
      LHS.submat(0, n, n-1, sizeLHS-1) = A.t();
    }
    mat RHS(sizeLHS, 1);
    // Initialising Nesterov-Todd scalings
    std::vector<std::map<std::string,mat> > WList;
    WList = initnts();
    // Initialising PDV for determining (first) interior point
    pdv->x = -q;
    pdv->y = b;
    pdv->z = cList.hmats;
    pdv = sxyz(pdv, LHS, RHS, WList);
    cps->set_pdv(*pdv);
  }
    */ 
  return cps;
}
