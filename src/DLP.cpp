#include "cccp.h"
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
  double term1 = 0.0, term2 = 0.0, ans;
  term1 = dot(b, pdv.y);
  term2 = sum(cList.sdot(pdv.z, cList.h));
  ans = -term1 - term2;
  return ans;
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
  Main routine for solving a Linear Program
*/
CPS* DLP::cps(CTRL& ctrl){
  // Initializing objects
  PDV* pdv = cList.initpdv(A.n_rows);
  PDV InitPrim, InitDual, KktSol;

  PDV* dpdv1 = cList.initpdv(A.n_rows);
  PDV* dpdv2 = cList.initpdv(A.n_rows);
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  bool checkRgap = false;
  int m = sum(cList.dims), n = cList.n, sizeLHS = A.n_rows + A.n_cols; 
  double resx, resx0, resy, resy0, resz, resz0, pres, dres, 
    ts, nrms, tz, nrmz, tt, tk, tm, pcost, dcost, gap, rgap = NA_REAL, 
    hresx, hresy, hresz, hz, by, cx, rt, pinfres, dinfres, 
    nomin, denom, dg = 1.0, dgi = 1.0, lg = 1.0, lgprd, dkdt = 0.0, mu, sigma, step;
  Rcpp::NumericVector state = cps->get_state();
  vec ss(5);
  mat rx, ry, rz, hrx, hry, hrz, Lambda, LambdaPrd, Ws3, Wh, Whz, dsdz;
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
  // Initialising LHS matrices
  LHS.zeros();
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(n, 0, sizeLHS-1, n-1) = A;
    LHS.submat(0, n, n-1, sizeLHS-1) = A.t();
  }
  // Computing initial values of PDV / CPS and scalings
  WList = cList.initnts();
  // Primal Start
  InitPrim.x = zeros(n, 1);
  InitPrim.y = b;
  InitPrim.z = cList.h;
  InitPrim = *(sxyz(&InitPrim, LHS, RHS, WList));
  InitPrim.s = -InitPrim.z;
  ts = cList.smss(InitPrim.s).max();
  // Dual Start
  InitDual.x = -q;
  InitDual.y = zeros(b.n_rows, 1);
  InitDual.z = zeros(cList.h.n_rows, 1);
  InitDual = *(sxyz(&InitDual, LHS, RHS, WList));
  tz = cList.smss(InitDual.z).max();
  // Initial Point
  pdv->x = InitPrim.x;
  pdv->y = InitDual.y;
  pdv->s = InitPrim.s;
  pdv->z = InitDual.z;
  nrms = sum(cList.snrm2(pdv->s));
  nrmz = sum(cList.snrm2(pdv->z));
  // Initial point optimal?
  gap = sum(cList.sdot(pdv->s, pdv->z));
  pcost = pobj(*pdv);
  dcost = -dot(b, pdv->y) - sum(cList.sdot(pdv->z, cList.h));
  if(pcost < 0.0) rgap = gap / (-pcost);
  if(dcost > 0.0) rgap = gap / dcost;
  if(!std::isnan(rgap)){
    checkRgap = (rgap <= rtol);
  } else {
    checkRgap = false;
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
    hresx = sqrt(dot(hrx, hrx));
    rx = hrx - q * pdv->tau;
    resx = sqrt(dot(rx, rx)) / pdv->tau;
    // Primal residuals
    hry = A * pdv->x;
    hresy = sqrt(dot(hry, hry));
    ry = hry - b * pdv->tau;
    resy = sqrt(dot(ry, ry)) / pdv->tau;
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
    if(!std::isnan(rgap)){
      checkRgap = (rgap <= rtol);
    } else {
      checkRgap = false;
    } 
    if((pres <= ftol) && (dres <= ftol) && ((gap <= atol) || checkRgap)){
      pdv->x = pdv->x / pdv->tau;
      pdv->y = pdv->y / pdv->tau;
      pdv->s = pdv->s / pdv->tau;
      pdv->z = pdv->z / pdv->tau;
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
      cps->set_pdv(*pdv);
      if(trace){
	Rcpp::Rcout << "Optimal solution found." << std::endl;
      }
      return cps;
    } else if((!std::isnan(pinfres)) && (pinfres <= ftol)){
      denom = -hz - by;
      pdv->x = mat(0,1);
      pdv->y = pdv->y / denom;
      pdv->s = mat(0,1);
      pdv->y = pdv->y / denom;
      tz = cList.smss(pdv->z).max();
      state["pobj"] = NA_REAL;
      state["dobj"] = 1.0;
      state["dgap"] = NA_REAL;
      state["rgap"] = NA_REAL;
      state["certp"] = pinfres;
      state["certd"] = NA_REAL;
      state["pslack"] = NA_REAL;
      state["dslack"] = -tz;
      cps->set_state(state);
      cps->set_status("primal infeasible");
      cps->set_niter(i);
      cps->set_pdv(*pdv); 
     if(trace){
	Rcpp::Rcout << "Certificate of primal infeasibility found." << std::endl;
      }
     return cps;
    } else if((!std::isnan(dinfres)) && (dinfres <= ftol)){
      denom = -cx;
      pdv->x = pdv->x / denom;
      pdv->y = mat(0,1);
      pdv->s = pdv->s / denom;
      pdv->z = mat(0,1);
      ts = cList.smss(pdv->s).max();
      state["pobj"] = -1.0;
      state["dobj"] = NA_REAL;
      state["dgap"] = NA_REAL;
      state["rgap"] = NA_REAL;
      state["certp"] = NA_REAL;
      state["certd"] = dinfres;
      state["pslack"] = -ts;
      state["dslack"] = NA_REAL;
      cps->set_state(state);
      cps->set_status("dual infeasible");
      cps->set_pdv(*pdv); 
      cps->set_niter(i);
      if(trace){
	Rcpp::Rcout << "Certificate of dual infeasibility found." << std::endl;
      }
     return cps;
    }
    // Computing initial scalings
    if(i == 0){
      WList = cList.ntsc(pdv->s, pdv->z);
      Lambda = cList.getLambda(WList);
      dg = sqrt(pdv->kappa / pdv->tau);
      dgi = sqrt(pdv->tau / pdv->kappa);
      lg = sqrt(pdv->tau * pdv->kappa);
    }
    LambdaPrd = cList.sprd(Lambda, Lambda);
    lgprd = lg * lg;
    // Solution step 1 (same for affine and combined solution)
    try{
      dpdv1->x = -q;
      dpdv1->y = b;
      dpdv1->z = cList.h;
      dpdv1 = sxyz(dpdv1, LHS, RHS, WList);
    } catch(std::runtime_error &ex){
      pdv->x = pdv->x / pdv->tau;
      pdv->y = pdv->y / pdv->tau;
      pdv->s = pdv->s / pdv->tau;
      pdv->z = pdv->z / pdv->tau;
      ts = cList.smss(pdv->s).max();
      tz = cList.smss(pdv->z).max();
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
      cps->set_status("unknown");
      cps->set_niter(i);
      cps->set_pdv(*pdv);
      if(trace){
	Rcpp::Rcout << "Terminated (singular KKT matrix)." << std::endl;
      }
      return cps;
    } catch(...){
      ::Rf_error("C++ exception (unknown reason)"); 
   }
    dpdv1->x = dpdv1->x * dgi;
    dpdv1->y = dpdv1->y * dgi;
    dpdv1->z = dpdv1->z * dgi;
    Wh = cList.ssnt(cList.h, WList, true, true);
    mu = dot(Lambda, Lambda) / (m + 1);
    sigma = 0.0;
    // Solving for affine and combined direction in two-round for-loop
    for(int ii = 0; ii < 2; ii++){
      dpdv2->s = LambdaPrd;
      dpdv2->kappa = lgprd;
      if(ii == 1){
	dpdv2->s = dpdv2->s + dsdz - OneE * sigma * mu;
	dpdv2->kappa = dpdv2->kappa + dkdt - sigma * mu;
      }
      dpdv2->x = (1.0 - sigma) * rx;
      dpdv2->y = (1.0 - sigma) * ry;
      dpdv2->z = (1.0 - sigma) * rz;
      dpdv2->tau = (1.0 - sigma) * rt;
      // Solving KKT-system
      dpdv2->s = cList.sinv(dpdv2->s, Lambda);
      dpdv2->s = -1.0 * dpdv2->s;
      Ws3 = cList.ssnt(dpdv2->s, WList, false, true);
      dpdv2->z = dpdv2->z + Ws3;
      dpdv2->z = -1.0 * dpdv2->z;
      dpdv2->y = -1.0 * dpdv2->y; 
      KktSol = *(sxyz(dpdv2, LHS, RHS, WList));

      // Combining solutions dpdv1 and dpdv2
      dpdv2->kappa = -1.0 * dpdv2->kappa / lg;
      dpdv2->tau = dpdv2->tau + dpdv2->kappa / dgi;
      cx = dot(q, KktSol.x);
      by = dot(b, KktSol.y);
      Whz = sum(cList.sdot(Wh, KktSol.z));
      nomin = (dgi * (dpdv2->tau + cx + by + Whz)).at(0,0);
      denom = 1.0 + sum(cList.sdot(dpdv1->z, dpdv1->z));
      dpdv2->tau = nomin / denom;
      dpdv2->x = KktSol.x + dpdv2->tau * dpdv1->x;
      dpdv2->y = KktSol.y + dpdv2->tau * dpdv1->y;
      dpdv2->z = KktSol.z + dpdv2->tau * dpdv1->z;
      dpdv2->s = dpdv2->s - dpdv2->z;
      dpdv2->kappa = dpdv2->kappa - dpdv2->tau;
      // ds o dz and dkappa * dtau for Mehrotra correction
      if(ii == 0){
	dsdz = cList.sprd(dpdv2->s, dpdv2->z);
	dkdt = dpdv2->kappa * dpdv2->tau; 
      }
      dpdv2->s = cList.sslb(dpdv2->s, Lambda, false);
      dpdv2->z = cList.sslb(dpdv2->z, Lambda, false); 
      ts = cList.smss(dpdv2->s).max();
      tz = cList.smss(dpdv2->z).max();
      tt = -dpdv2->tau / lg;
      tk = -dpdv2->kappa / lg;
      ss << 0.0 << ts << tz << tt << tk << endr;
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
	sigma = pow((1.0 - step), 3.0);
      }
    } // end ii-loop

    // Updating x, y; s and z (in current scaling)
    pdv->x = pdv->x + step * dpdv2->x;
    pdv->y = pdv->y + step * dpdv2->y;

    dpdv2->s = cList.SorZupdate(dpdv2->s, Lambda, step);
    dpdv2->z = cList.SorZupdate(dpdv2->z, Lambda, step);

    // Updating NT-scaling and Lagrange Multipliers
    WList = cList.ntsu(dpdv2->s, dpdv2->z, WList);
    Lambda = cList.getLambda(WList);
    dg = dg * sqrt(1.0 - step * tk) / sqrt(1.0 - step * tt);
    dgi = 1 / dg;
    lg = lg * sqrt(1 - step * tt) * sqrt(1 - step * tk);
    pdv->s = cList.ssnt(Lambda, WList, false, true);
    pdv->z = cList.ssnt(Lambda, WList, true, false);
    pdv->kappa = lg / dgi;
    pdv->tau = lg * dgi;
    gap = pow(sqrt(dot(Lambda, Lambda)) / pdv->tau, 2.0);
  } // end i-loop

  pdv->x = pdv->x / pdv->tau;
  pdv->y = pdv->y / pdv->tau;
  pdv->s = pdv->s / pdv->tau;
  pdv->z = pdv->z / pdv->tau;
  ts = cList.smss(pdv->s).max();
  tz = cList.smss(pdv->z).max();
  cps->set_pdv(*pdv);
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
