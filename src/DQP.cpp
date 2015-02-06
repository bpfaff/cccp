#include "cccp.h"
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
  if(A.n_rows > 0){
    term1 = pdv.y.t() * (A * pdv.x - b);
  }
  // dobj term for inequality constraints
  if(cList.K > 0){
    term2 = dot(pdv.z, cList.G * pdv.x - cList.h);
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
Centrality Residuals
*/
mat DQP::rcent(PDV& pdv){
  mat ans(cList.G.n_rows, 1);

  ans = pdv.s + cList.G * pdv.x - cList.h;

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
    Gz = cList.G.t() * pdv.z;
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
    nomin = cList.snrm2(rz);
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

  nomin = norm(rdual(pdv));
  denom = std::max(1.0, norm(q));
  ans = nomin / denom;

  return ans;
}
/*
  Solving 'KKT-System'
*/
PDV* DQP::sxyz(PDV* pdv, mat LHS, mat RHS, std::vector<std::map<std::string,mat> > WList){
  int n = P.n_cols;
  mat lhs1, rhs1, ans;

  lhs1 = cList.gwwg(WList);
  LHS.submat(0, 0, n-1, n-1) = P + lhs1;
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
CPS* DQP::cps(CTRL& ctrl){
  // Initialising object
  PDV* pdv = cList.initpdv(A.n_rows);
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  Rcpp::List params(ctrl.get_params());
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
    double ftol = Rcpp::as<double>(params["feastol"]);
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
  // Defining variables used in iterations
  bool trace = Rcpp::as<bool>(params["trace"]);
  bool checkRgap = false;
  int m = sum(cList.dims), n = P.n_cols, sizeLHS = A.n_rows + A.n_cols; 
  int maxiters = Rcpp::as<int>(params["maxiters"]);
  double resx, resx0, resy, resy0, resz, resz0, 
    pcost, dcost, gap, rgap = NA_REAL, pres, dres, 
    ts, nrms, tz, nrmz, tm, step, mu, sigma, dsdz;
  double atol = Rcpp::as<double>(params["abstol"]);
  double ftol = Rcpp::as<double>(params["feastol"]);
  double rtol = Rcpp::as<double>(params["reltol"]);
  double sadj = Rcpp::as<double>(params["stepadj"]);
  vec ss(3);
  mat LHS(sizeLHS, sizeLHS);
  mat RHS(sizeLHS, 1);
  mat rx, ry, rz, Lambda, LambdaPrd, Ws3;
  mat OneE = cList.sone();
  std::vector<std::map<std::string,mat> > WList;
  PDV* dpdv = cList.initpdv(A.n_rows);
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
  // Initialising Nesterov-Todd scalings
  WList = cList.initnts();
  // Initialising PDV for determining (first) interior point
  pdv->x = -q;
  pdv->y = b;
  pdv->z = cList.h;
  pdv = sxyz(pdv, LHS, RHS, WList);
  pdv->s = -1.0 * (pdv->z);
  ts = cList.smss(pdv->s).max();
  nrms = sum(cList.snrm2(pdv->s));
  tz = cList.smss(pdv->z).max();
  nrmz = sum(cList.snrm2(pdv->z));
  if(ts >= -1e-8 * std::max(1.0, nrms)){
    pdv->s = cList.sams1(pdv->s, ts);
  }
  if(tz >= -1e-8 * std::max(1.0, nrmz)){
    pdv->z = cList.sams1(pdv->z, tz);
  }
  // Duality gap for initial solution
  gap = sum(cList.sdot(pdv->s, pdv->z));
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
    resz = cList.snrm2(rz);
    // Statistics for stopping criteria
    pcost = pobj(*pdv);
    dcost = pcost + dot(pdv->y, ry) + cList.sdot(pdv->z, rz).at(0, 0) - gap;
    rgap = NA_REAL;
    if(pcost < 0.0) rgap = gap / (-pcost);
    if(dcost > 0.0) rgap = gap / dcost;
    pres = std::max(resy / resy0, resz / resz0); 
    dres = resx / resx0;
    // Tracing status quo of IPM
    if(trace){
      Rcpp::Rcout << "Iteration: " << i << std::endl;
      Rcpp::Rcout << "pobj: " << pcost << std::endl;
      Rcpp::Rcout << "dobj: " << dcost << std::endl;
      Rcpp::Rcout << "pinf: " << pres << std::endl;
      Rcpp::Rcout << "dinf: " << dres << std::endl;
      Rcpp::Rcout << "dgap: " << gap << std::endl;
      Rcpp::Rcout << std::endl;
    }
    // Checking convergence
    if(!std::isnan(rgap)){
      checkRgap = (rgap <= rtol);
    } else {
      checkRgap = false;
    } 
    if((pres <= ftol) && (dres <= ftol) && ((gap <= atol) || checkRgap)){
      cps->set_pdv(*pdv);
      cps->set_sidx(cList.sidx);
      ts = cList.smss(pdv->s).max();
      tz = cList.smss(pdv->z).max();
      state["pobj"] = pobj(*pdv);
      state["dobj"] = dobj(*pdv);
      state["dgap"] = gap;
      state["certp"] = certp(*pdv);
      state["certd"] = certd(*pdv);
      state["pslack"] = -ts;
      state["dslack"] = -tz;
      if(!std::isnan(rgap)){
	state["rgap"] = rgap;
      }

      cps->set_state(state);
      cps->set_status("optimal");
      cps->set_niter(i);
      if(trace){
	Rcpp::Rcout << "Optimal solution found." << std::endl;
      }
      return cps;
    }
    // Computing initial scalings
    if(i == 0){
      WList = cList.ntsc(pdv->s, pdv->z);
    }
    Lambda = cList.getLambda(WList);
    LambdaPrd = cList.sprd(Lambda, Lambda);
    mu = gap / m;
    sigma = 0.0;
    // Solving for affine and combined direction in two-round for-loop
    for(int ii = 0; ii < 2; ii++){
      dpdv->x = -rx;
      dpdv->y = -ry;
      dpdv->z = -rz;
      // Solving KKT-System
      dpdv->s = -LambdaPrd + OneE * sigma * mu;
      dpdv->s = cList.sinv(dpdv->s, Lambda);
      Ws3 = cList.ssnt(dpdv->s, WList, false, true);
      dpdv->z = dpdv->z - Ws3;
      dpdv = sxyz(dpdv, LHS, RHS, WList);
      dpdv->s = dpdv->s - dpdv->z;
      // ds o dz for Mehrotra correction
      dsdz = sum(cList.sdot(dpdv->s, dpdv->z));
      dpdv->s = cList.sslb(dpdv->s, Lambda, false);
      dpdv->z = cList.sslb(dpdv->z, Lambda, false);

      ts = cList.smss(dpdv->s).max();
      tz = cList.smss(dpdv->z).max();
      ss << 0.0 << ts << tz << endr;
      tm = ss.max();
      if(tm == 0.0){
	step = 1.0;
      } else {
	if(ii == 0) {
	  step = std::min(1.0, 1.0 / tm);
	} else {
	  step = std::min(1.0, sadj / tm);
	}
      }
      if(ii == 0){
	sigma = pow(std::min(1.0, std::max(0.0, 1.0 - step + dsdz / gap * std::pow(step, 2.0))), 3.0);
      }
    } // end ii-loop

    // Updating x, y; s and z (in current scaling)
    pdv->x = pdv->x + step * dpdv->x;
    pdv->y = pdv->y + step * dpdv->y;

    dpdv->s = cList.SorZupdate(dpdv->s, Lambda, step);
    dpdv->z = cList.SorZupdate(dpdv->z, Lambda, step);

    // Updating NT-scaling and Lagrange Multipliers
    WList = cList.ntsu(dpdv->s, dpdv->z, WList);
    Lambda = cList.getLambda(WList);
    pdv->s = cList.ssnt(Lambda, WList, false, true);
    pdv->z = cList.ssnt(Lambda, WList, true, false);
    gap = sum(cList.sdot(Lambda, Lambda));
  } // end i-loop

  // Preparing result for non-convergence in maxiters iterations
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  state["pobj"] = pobj(*pdv);
  state["dobj"] = dobj(*pdv);
  state["dgap"] = gap;
  state["certp"] = certp(*pdv);
  state["certd"] = certd(*pdv);
  state["pslack"] = -ts;
  state["dslack"] = -tz;
  if(!std::isnan(rgap)){
    state["rgap"] = rgap;
  }
  cps->set_state(state);
  cps->set_niter(maxiters);

  if((state["certp"] <= ftol) && (state["certd"] <= ftol)){
    cps->set_status("optimal");
  } else {
    if(trace){
      Rcpp::Rcout << "Optimal solution not determined in " << maxiters << " iteration(s)." << std::endl;
    }
    cps->set_status("unknown");
  }
  return cps;
}
