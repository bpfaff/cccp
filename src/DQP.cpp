#include "cccp3.h"
/*
 *
 * Methods for Quadratic Programs
 *
*/
using namespace arma;
const int TraceFieldWidth = 10;
const int TracePrintPrecs = 5;
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
      term2 = term2 + dot(pdv.z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), \
	(cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * pdv.x - \
	 cList.h(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all)));
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
  mat ans(cList.G.n_rows, 1);

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
      W.insert(std::pair<std::string,mat>("dnl", ans));
      W.insert(std::pair<std::string,mat>("dnli", ans));
      ans = zeros(cList.dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cList.cone[i] == "NNOC"){
      ans = ones(cList.dims[i],1);
      W.insert(std::pair<std::string,mat>("d", ans));
      W.insert(std::pair<std::string,mat>("di", ans));
      ans = zeros(cList.dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cList.cone[i] == "SOCC"){
      ans = ones(1,1);
      W.insert(std::pair<std::string,mat>("beta", ans));
      ans = zeros(cList.dims[i],1);
      ans.at(0,0) = 1.0;
      W.insert(std::pair<std::string,mat>("v", ans));
      ans = zeros(cList.dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cList.cone[i] == "PSDC"){
      ans = eye(cList.dims[i],cList.dims[i]);
      ans.reshape(cList.dims[i] * cList.dims[i], 1);
      W.insert(std::pair<std::string,mat>("r", ans));
      W.insert(std::pair<std::string,mat>("rti", ans));
      ans = zeros(cList.dims[i] * cList.dims[i], 1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    }
    WList.push_back(W);
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
  // Computing fixed values
  //int m = sum(cList.dims);
  double resx0 = std::max(1.0, norm(q));
  double resy0 = std::max(1.0, norm(b));
  double resz0 = std::max(1.0, snrm2(cList.h));
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
  pdv->z = cList.h;
  pdv = sxyz(pdv, LHS, RHS, WList);
  pdv->s = -1.0 * (pdv->z);
  double ts = smss(pdv->s).max();
  double nrms = sum(snrm2(pdv->s));
  double tz = smss(pdv->z).max();
  double nrmz = sum(snrm2(pdv->z));
  if(ts >= -1e-8 * std::max(1.0, nrms)){
    pdv->s = sams1(pdv->s, ts);
  }
  if(tz >= -1e-8 * std::max(1.0, nrmz)){
    pdv->z = sams1(pdv->z, tz);
  }
  // Defining variables used in iterations
  bool trace = ctrl.get_trace(), checkRgap = false;
  int maxiters = ctrl.get_maxiters();
  double resx, resy, resz, pcost, dcost, gap, rgap = NA_REAL, pres, dres;
  double atol = ctrl.get_abstol();
  double ftol = ctrl.get_feastol();
  double rtol = ctrl.get_reltol();
  mat rx, ry, rz;
  // Duality gap for initial solution
  gap = sum(sdot(pdv->s, pdv->z));
  cps->set_pdv(*pdv);
  //
  // Starting iterations
  //
  for(int i = 0; i < maxiters; i++){
    // Dual Residuals
    rx = rdual(*pdv);
    resx = norm(rx);
    // Primal Residuals
    ry = rprim(*pdv);
    resy = norm(ry);
    // Central Residuals 
    rz = rcent(*pdv);
    resz = snrm2(rz);
    // Statistics for stopping criteria
    pcost = pobj(*pdv);
    dcost = pcost + dot(pdv->y, ry) + sdot(pdv->z, rz).at(0, 0) - gap;
    if(pcost < 0.0) rgap = gap / (-pcost);
    if(dcost > 0.0) rgap = gap / dcost;
    pres = std::max(resy / resy0, resz / resz0); 
    dres = resx / resx0;
    // Tracing status quo of IPM
    if(trace){
      Rcpp::Rcout << "Iteration: " << i + 1 << std::endl;
      Rcpp::Rcout << std::setiosflags(ios::left) 
		  << std::setw(TraceFieldWidth) << "pobj" 
		  << std::setw(TraceFieldWidth) << "dobj" 
		  << std::setw(TraceFieldWidth) << "pinf" 
		  << std::setw(TraceFieldWidth) << "dinf" 
		  << std::setw(TraceFieldWidth) << "dgap" 
		  << std::endl;
      Rcpp::Rcout << std::setiosflags(ios::left) 
		  << std::setw(TraceFieldWidth) << std::setprecision(TracePrintPrecs) << pcost 
		  << std::setw(TraceFieldWidth) << std::setprecision(TracePrintPrecs) << dcost 
		  << std::setw(TraceFieldWidth) << std::setprecision(TracePrintPrecs) << pres 
		  << std::setw(TraceFieldWidth) << std::setprecision(TracePrintPrecs) << dres 
		  << std::setw(TraceFieldWidth) << std::setprecision(TracePrintPrecs) << gap 
		  << std::endl;
    }
    // Checking convergence
    if(!std::isnan(rgap)){
      checkRgap = (rgap <= rtol);
    }
    if((pres <= ftol) && (dres <= ftol) && ((gap <= atol) || checkRgap)){
      cps->set_pdv(*pdv);

      ts = smss(pdv->s).max();
      tz = smss(pdv->z).max();
      state["pobj"] = pobj(*pdv);
      state["dobj"] = dobj(*pdv);
      state["dgap"] = gap;
      state["certp"] = certp(*pdv);
      state["certd"] = certd(*pdv);
      state["pslack"] = -ts;
      state["dslack"] = -tz;
      cps->set_state(state);
      cps->set_status("optimal");
      cps->set_niter(i + 1);
      if(trace){
	Rcpp::Rcout << "Optimal solution found." << std::endl;
      }

    }
    /*
      Computing initial scalings
    */


  } 


  // Preparing result (not complete, yet)
  state["pobj"] = pobj(*pdv);
  state["dobj"] = dobj(*pdv);
  state["dgap"] = gap;
  state["certp"] = certp(*pdv);
  state["certd"] = certd(*pdv);
  state["pslack"] = -ts;
  state["dslack"] = -tz;
  cps->set_state(state);
  cps->set_niter(maxiters);

  if((state["certp"] <= ftol) && (state["certd"] <= ftol)){
    cps->set_status("optimal");
  } else {
    cps->set_status("unknown");
  }
  if(trace){
    Rcpp::Rcout << "Optimal solution not determined in " << maxiters << " iteration(s)." << std::endl;
  }

  return cps;
}
