#include "cccp3.h"
/*
 *
 * Methods for Linear Programs
 *
*/
using namespace arma;
/*
Primal objective
*/
double DLP::pobj(PDV& pdv){
  double ans = dot(pdv.x, q);
  return ans;
}
/*
Dual objective
*/
double DLP::dobj(PDV& pdv){
  double term1 = 0.0, term2 = 0.0;
  term1 = dot(b, pdv.y);
  term2 = sum(cList.sdot(pdv.z, cList.h));
  return -term1 - term2;
}
/*
Primal Residuals
*/
mat DLP::rprim(PDV& pdv){
  int p = A.n_rows;
  mat ans(p,1);
  ans.zeros();

  ans = b - A * pdv.x;

  return ans;
}
/*
Centrality Residuals
*/
mat DLP::rcent(PDV& pdv){
  mat ans(cList.G.n_rows, 1);

  ans = pdv.s + cList.G * pdv.x - cList.h;

  return ans;
}
/*
Dual Residuals
*/
mat DLP::rdual(PDV& pdv){
  int n = q.n_rows;
  mat Gz(n,1);
  mat Ay(n,1);
  mat ans(n,1);
  Gz.zeros();
  Ay.zeros();
  ans.zeros();

  if(cList.K > 0){
    Gz = cList.G.t() * pdv.z;
  }
  if(A.n_rows > 0){
    Ay = A.t() * pdv.y;
  }
  ans = q + Gz + Ay;

  return ans;
}
/*
Certificate of primal infeasibilty
*/
double DLP::certp(PDV& pdv){
  double nomin, denom, ans1 = 0.0, ans2 = 0.0, ans = 0.0;

  nomin = norm(rprim(pdv));
  denom = std::max(1.0, norm(b));
  ans1 = nomin / denom;

  if(cList.K > 0){
    mat rz;
    rz = rcent(pdv);
    nomin = 0.0;
    denom = std::max(1.0, norm(q));
    nomin = cList.snrm2(rz);
    ans2 = nomin / denom;
  }
  ans = std::max(ans1, ans2);

  return ans;
}
/*
Certificate of dual infeasibilty
*/
double DLP::certd(PDV& pdv){
  double nomin, denom, ans;

  nomin = norm(rdual(pdv));
  denom = std::max(1.0, norm(q));
  ans = nomin / denom;

  return ans;
}
/*
  Solving 'KKT-System'
*/
PDV* DLP::sxyz(PDV* pdv, mat LHS, mat RHS, std::vector<std::map<std::string,mat> > WList){
  int n = q.n_rows;
  mat lhs1, rhs1, ans;

  lhs1 = cList.gwwg(WList);
  LHS.submat(0, 0, n-1, n-1) = lhs1;
  rhs1 = cList.gwwz(WList, pdv->z);
  RHS.submat(0, 0, n - 1, 0) = pdv->x + rhs1;
  if(pdv->y.n_rows > 0){
    RHS.submat(n, 0, RHS.n_rows - 1, 0) = pdv->y;
  }
  ans = solve(LHS, RHS);
  pdv->x = ans.submat(0, 0, n - 1, 0);
  if(pdv->y.n_rows > 0){
    pdv->y = ans.submat(n, 0, RHS.n_rows - 1, 0);
  }
  pdv->z = cList.G * pdv->x - pdv->z;
  pdv->z = cList.ssnt(pdv->z, WList, true, true);

  return pdv;
}
/*
  Main routine for solving a Quadratic Program
*/
CPS* DLP::cps(CTRL& ctrl){
  // Initializing objects
  PDV* pdv = cList.initpdv(A.n_rows);
  PDV* InitPrim = cList.initpdv(A.n_rows);
  PDV* InitDual = cList.initpdv(A.n_rows);
  PDV* dpdv1 = cList.initpdv(A.n_rows);
  PDV* dpdv2 = cList.initpdv(A.n_rows);
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  bool checkRgap = false;
  int m = sum(cList.dims), n = cList.n, sizeLHS = A.n_rows + A.n_cols; 
  double resx, resx0, resy, resy0, resz, resz0, pres, dres, 
    ts, nrms, tz, nrmz, pcost, dcost, gap, rgap, hresx, hresy, hresz,
    hz, by, cx, rt, pinfres, dinfres;
  Rcpp::NumericVector state = cps->get_state();
  vec ss(3), sdv(3), eval;
  mat rx, ry, rz, hrx, hry, hrz;
  mat OneE = cList.sone();
  mat LHS(sizeLHS, sizeLHS);
  mat RHS(sizeLHS, 1);
  std::vector<std::map<std::string,mat> > WList;
  // Setting control parameters
  Rcpp::List params(ctrl.get_params());
  bool trace = Rcpp::as<bool>(params["trace"]);
  int maxiters = Rcpp::as<int>(params["maxiters"]);
  double atol = Rcpp::as<double>(params["abstol"]);
  double ftol = Rcpp::as<double>(params["feastol"]);
  double rtol = Rcpp::as<double>(params["reltol"]);
  double sadj = Rcpp::as<double>(params["stepadj"]);
  // Computing fixed values
  resx0 = std::max(1.0, norm(q));
  resy0 = std::max(1.0, norm(b));
  resz0 = std::max(1.0, cList.snrm2(cList.h));
  // Initialising LHS and RHS matrices
  LHS.zeros();
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(n, 0, sizeLHS-1, n-1) = A;
    LHS.submat(0, n, n-1, sizeLHS-1) = A.t();
  }
  // Computing initial values of PDV / CPS and scalings
  WList = cList.initnts();
  // Primal Start
  pdv->x = pdv->x.zeros();
  pdv->y = b;
  pdv->z = cList.h;
  InitPrim = sxyz(pdv, LHS, RHS, WList);
  InitPrim->s = -InitPrim->z;
  ts = cList.smss(InitPrim->s).max();
  // Dual Start
  pdv->x = -q;
  pdv->y = b.zeros();
  pdv->z = cList.h.zeros();
  InitDual = sxyz(pdv, LHS, RHS, WList);
  tz = cList.smss(InitDual->z).max();
  // Initial Point
  pdv->x = InitPrim->x;
  pdv->y = InitDual->y;
  pdv->s = InitPrim->s;
  pdv->z = InitDual->z;
  nrms = sum(cList.snrm2(pdv->s));
  nrmz = sum(cList.snrm2(pdv->z));
  // Initial point optimal?
  gap = sum(cList.sdot(pdv->s, pdv->z));
  pcost = pobj(*pdv);
  dcost = -dot(b, pdv->y) - cList.sdot(pdv->z, cList.h).at(0, 0);
  if(pcost < 0.0) rgap = gap / (-pcost);
  if(dcost > 0.0) rgap = gap / dcost;
  if(!std::isnan(rgap)){
    checkRgap = (rgap <= rtol);
  }
  if((ts <= ftol) && (tz <= ftol) && ((gap <= atol) || checkRgap)){
    // Dual Residuals
    rx = rdual(*pdv);
    resx = norm(rx);
    // Primal Residuals
    ry = rprim(*pdv);
    resy = norm(ry);
    // Central Residuals 
    rz = rcent(*pdv);
    resz = cList.snrm2(rz);
    pres = std::max(resy / resy0, resz / resz0); 
    dres = resx / resx0;
    cps->set_pdv(*pdv);
    cps->set_sidx(cList.sidx);
    state["pobj"] = pobj(*pdv);
    state["dobj"] = dobj(*pdv);
    state["dgap"] = gap;
    state["certp"] = pres;
    state["certd"] = dres;
    state["pslack"] = -ts;
    state["dslack"] = -tz;
    if(!std::isnan(rgap)){
      state["rgap"] = rgap;
    }
    cps->set_state(state);
    cps->set_status("optimal");
    cps->set_niter(0);
    if(trace){
      Rcpp::Rcout << "Optimal solution found." << std::endl;
    }
    return cps;
  } 
  if(ts >= -1e-8 * std::max(1.0, nrms)){
    pdv->s = cList.sams1(pdv->s, ts);
  }
  if(tz >= -1e-8 * std::max(1.0, nrmz)){
    pdv->z = cList.sams1(pdv->z, tz);
  }
  // Duality gap for initial solution
  gap = sum(cList.sdot(pdv->s, pdv->z));
  //
  // Starting iterations
  //
  for(int i = 0; i < maxiters; i++){
    // Evaluate residuals, gap and stopping criteria
    // Dual residuals
    hrx = -(A.t() * pdv->y) - cList.G.t() * pdv->z;
    hresx = norm(hrx);
    rx = hrx - q * pdv->tau;
    resx = norm(rx) / pdv->tau;
    // Primal residuals
    hry = A * pdv->x;
    hresy = norm(hry);
    ry = hry - b * pdv->tau;
    resy = norm(ry) / pdv->tau;
    // Centrality residuals
    hrz = pdv->s + cList.G * pdv->x;
    hresz = cList.snrm2(hrz);
    rz = hrz - cList.h * pdv->tau;
    resz = cList.snrm2(rz) / pdv->tau;
    // Self-dual residuals
    hz = sum(cList.sdot(cList.h, pdv->z));
    by = dot(b, pdv->y);
    cx = dot(q, pdv->x);
    rt = pdv->kappa + cx + by + hz;
    // Statistics for stopping criteria
    pcost = cx / pdv->tau;
    dcost = -(by + hz) / pdv->tau;
    rgap = NA_REAL;
    if(pcost < 0.0) rgap = gap / (-pcost);
    if(dcost > 0.0) rgap = gap / dcost;
    pres = std::max(resy / resy0, resz / resz0); 
    dres = resx / resx0;
    if(hz + by < 0.0){
      pinfres = hresx / resx0 / (-hz - by);
    } else {
      pinfres = NA_REAL;
    }
    if(cx < 0.0){
      dinfres = std::max(hresy / resy0, hresz / resz0) / (-cx);
    } else {
      dinfres = NA_REAL;
    }
    if(trace){
      Rcpp::Rcout << "Iteration: " << i << std::endl;
      Rcpp::Rcout << "pobj: " << pcost << std::endl;
      Rcpp::Rcout << "dobj: " << dcost << std::endl;
      Rcpp::Rcout << "pinf: " << pres << std::endl;
      Rcpp::Rcout << "dinf: " << dres << std::endl;
      Rcpp::Rcout << "dgap: " << gap << std::endl;
      Rcpp::Rcout << "k/t : " << pdv->kappa / pdv->tau << std::endl;
      Rcpp::Rcout << std::endl;
    }
    // Checking convergence / infeasibilities


  }

  cps->set_pdv(*pdv);

  return cps;
}
