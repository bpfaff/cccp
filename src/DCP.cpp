#include "cccp.h"
/*
 *
 * Methods for Convex Programs with nonlinear constraints
 *
*/
using namespace arma;
/*
Primal objective
*/
double DCP::pobj(PDV& pdv){
  double ans = pdv.x.at(pdv.x.n_rows - 1, 0);
  return ans;
}
/*
Dual objective
*/
double DCP::dobj(PDV& pdv){
  double term1 = 0.0, term2 = 0.0, term3 = 0.0, ans;
  term1 = pdv.x(pdv.x.n_rows - 1, 0);
  if(cList.K > 1){
    for(int i = 1; i < cList.K; i++){
      term2 += dot(pdv.z(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all), 
		   cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * 
		   pdv.x - 
		   cList.h(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all));
    }
  }
  term2 += dot(pdv.z(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all), 
	       cList.h(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all));
  term3 = dot(pdv.y, (A * pdv.x - b));
  ans = term1 + term2 + term3;

  return ans;
}
/*
Primal Residuals
*/
mat DCP::rprim(PDV& pdv){
  int p = A.n_rows;
  mat ans(p, 1);
  ans.zeros();

  ans = b - A * pdv.x;

  return ans;
}
/*
Centrality Residuals
*/
mat DCP::rcent(PDV& pdv){
  mat ans(cList.G.n_rows, 1);

  ans(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all) = 
    pdv.s(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all) + 
    cList.h(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all);
 
 if(cList.K > 1){
    for(int i = 1; i < cList.K; i++){
      ans(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) = 
	pdv.s(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) +
	cList.G(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all) * pdv.x -
	cList.h(span(cList.sidx.at(i, 0), cList.sidx.at(i, 1)), span::all);
    }
  }

  return ans;
}
/*
Dual Residuals
*/
mat DCP::rdual(PDV& pdv){
  int n = x0.n_rows;
  mat Gz(n,1);
  mat Ay(n,1);
  mat ans(n,1);
  Gz.zeros();
  Ay.zeros();
  ans.zeros();

  Gz = cList.G.t() * pdv.z;
  
  if(A.n_rows > 0){
    Ay = A.t() * pdv.y;
  }
  ans = Gz + Ay;
  ans.at(ans.n_rows - 1, 0) += 1.0;

  return ans;
}
/*
Certificate of primal infeasibilty
*/
double DCP::certp(PDV& pdv){
  double nomin, denom, ans1 = 0.0, ans2 = 0.0, ans = 0.0;
  mat rz;

  nomin = norm(rprim(pdv));
  denom = std::max(1.0, norm(b));
  ans1 = nomin / denom;

  rz = rcent(pdv);
  ans2 = cList.snrm2(rz);

  ans = std::max(ans1, ans2);

  return ans;
}
/*
Certificate of dual infeasibilty
*/
double DCP::certd(PDV& pdv){
  double ans;

  ans = norm(rdual(pdv));

  return ans;
}
/*
  Solving 'KKT-System'
*/
PDV* DCP::sxyz(PDV* pdv, mat LHS, std::vector<std::map<std::string,mat> > WList){

  int n = x0.n_rows;
  int ne = n - 1;
  int mnl = cList.dims(0);
  int K = cList.K;
  double a = pdv->z.at(0, 0); // Slack with respect to f0
  double x1 = pdv->x.at(pdv->x.n_rows, 0); // Epigraph-variable 't'
  mat ux = mat(ne, 1), uz = pdv->z, RHS(ne + A.n_rows, 1), ans(ne + A.n_rows, 1);
  CONEC cEpi = cList;
  std::vector<std::map<std::string,mat> > WEpi = WList;
  
  cEpi.n -= 1;
  cEpi.G.set_size(cList.G.n_rows, ne);
  cEpi.G = cList.G(span::all, span(0, ne - 1)); // removing last column
  ux = pdv->x(span(0, ne - 1), span::all);
  ux = ux + pdv->x.at(n - 1, 0) * cEpi.G.row(0).t();

  // Distinguishing four cases:
  // mnl == 1 and K == 1: only f0 and no cone constraints
  // mnl == 1 and K > 1 : only f0 and cone constraints
  // mnl > 1 and K == 1 : f0 and other nonlinear constraints; no cone constraints
  // mnl > 1 and K > 1 : f0 and other nonlinear constraints and cone constraints
  // Problem to be solved is reduced to x0

  // mnl == 1 and K == 1
  if((mnl == 1) && (K == 1)){
    // upper LHS is only Hessian
    RHS.submat(0, 0, ne - 1, 0) = ux;
  }
  // mnl == 1 and K > 1
  if((mnl == 1) && (K > 1)){
    WEpi.erase(WEpi.begin());
    cEpi.K -= 1;
    cEpi.G = cEpi.G(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    cEpi.h = cEpi.h(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    cEpi.sidx = cEpi.sidx(span(1, cEpi.sidx.n_rows - 1), span::all);
    cEpi.sidx -= 1;
    cEpi.sidx.at(0, 0) = 0;
    cEpi.cone.erase(cEpi.cone.begin());
    cEpi.dims.shed_row(0);
    LHS.submat(0, 0, ne - 1, ne - 1) += cEpi.gwwg(WEpi);
    uz = uz(span(1, uz.n_rows - 1), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  // mnl > 1 and K == 1
  if((mnl > 1) && (K == 1)){
    WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
    WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
    cEpi.dims(0) -= 1;
    cEpi.sidx(0, 1) -= 1;
    cEpi.G = cEpi.G(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    cEpi.h = cEpi.h(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    LHS.submat(0, 0, ne - 1, ne - 1) += cEpi.gwwg(WEpi);
    uz = uz(span(1, uz.n_rows - 1), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  // mnl > 1 and K > 1
  if((mnl > 1) && (K > 1)){
    WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
    WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
    cEpi.dims(0) -= 1;
    cEpi.sidx = cEpi.sidx - 1;
    cEpi.sidx(0, 0) = 0;
    cEpi.G = cEpi.G(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    cEpi.h = cEpi.h(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
    LHS.submat(0, 0, ne - 1, ne - 1) += cEpi.gwwg(WEpi);
    uz = uz(span(1, uz.n_rows - 1), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  if(pdv->y.n_rows > 0){
    RHS.submat(ne, 0, RHS.n_rows - 1, 0) = pdv->y;
  }
  // Solving KKT-system
  ans = solve(LHS, RHS);
  // Preparing pdv
  pdv->x.submat(0, 0, ne - 1, 0) = ans.submat(0, 0, ne - 1, 0);
  if(pdv->y.n_rows > 0){
    pdv->y = ans.submat(ne, 0, RHS.n_rows - 1, 0);
  }
  uz = cEpi.G * pdv->x.submat(0, 0, ne - 1, 0) - uz;
  if((mnl > 1) || (K > 1)){
    pdv->z(span(1, pdv->z.n_rows - 1), span::all) = cEpi.ssnt(uz, WEpi, true, true);
  }
  pdv->z.at(0, 0) = -pdv->x.at(pdv->x.n_rows - 1, 0) * WList[0]["dnl"].at(0, 0);
  x1 = dot(cList.G.submat(0, 0, 0, ne - 1), pdv->x.submat(0, 0, ne - 1, 0)) + 
    pow(WList[0]["dnl"].at(0, 0), 2) * pdv->x.at(n - 1, 0) - a;
  pdv->x.at(n - 1, 0) = x1;
  return pdv;
}
/*
  Main routine for solving a Convex Program with nonlinear constraints
*/
CPS* DCP::cps(CTRL& ctrl){
  // Initializing objects
  PDV* pdv = cList.initpdv(A.n_rows);
  PDV* dpdv = cList.initpdv(A.n_rows);
  pdv->x = x0;
  Rcpp::List nF(nList[0]);
  Rcpp::List gF(nList[1]);
  Rcpp::List hF(nList[2]);

  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  Rcpp::NumericVector state = cps->get_state();
  bool checkRgap = false, backTrack;
  int m = sum(cList.dims), mnl = cList.dims(0), n = cList.n, 
    ne = n - 1, sizeLHS = A.n_rows + A.n_cols - 1;
  double gap = m, resx, resy, resz, pcost, dcost, rgap = NA_REAL, 
    pres, dres, pres0 = 1.0, dres0 = 1.0, sigma, mu, ts, tz, tm, step;
  vec ss(3), Fval(mnl);
  mat H = zeros(ne, ne), rx, ry, rz, Lambda, LambdaPrd, Ws3, x;
  mat OneE = cList.sone();
  mat LHS(sizeLHS, sizeLHS);
  // Initialising LHS matrices
  LHS.zeros();
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(ne, 0, sizeLHS - 1, ne - 1) = A(span::all, span(0, ne - 1));
    LHS.submat(0, ne, ne - 1, sizeLHS - 1) = A(span::all, span(0, ne - 1)).t();
  }
  std::vector<std::map<std::string,mat> > WList;
  // Setting control parameters
  Rcpp::List params(ctrl.get_params());
  bool trace = Rcpp::as<bool>(params["trace"]);
  int maxiters = Rcpp::as<int>(params["maxiters"]);
  double atol = Rcpp::as<double>(params["abstol"]);
  double ftol = Rcpp::as<double>(params["feastol"]);
  double rtol = Rcpp::as<double>(params["reltol"]);
  double sadj = Rcpp::as<double>(params["stepadj"]);
  double beta = Rcpp::as<double>(params["beta"]);
  //
  // Starting iterations
  //
  for(int i = 0; i < maxiters; i++){
    H.zeros();
    for(int j = 0; j < mnl; j++){
      // Setting f to first mnl-rows of h-matrix
      cList.h(j, 0) = feval(pdv->x(span(0, ne - 1), span::all), nF[j]);
      // Setting Df to first mnl-rows of G-matrix
      cList.G(j, span(0, ne - 1)) = geval(pdv->x(span(0, ne - 1), span::all), gF[j]).t();
      // Computing Hessian
      H += pdv->z.at(j, 0) * heval(pdv->x(span(0, ne - 1), span::all), hF[j]);
    }
    cList.h(0, 0) = cList.h(0, 0) - pdv->x.at(n - 1, 0);
    // Computing gap
    gap = sum(cList.sdot(pdv->s, pdv->z));
    // Computing residuals
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
    dcost = pcost + dot(ry, pdv->y) + sum(cList.sdot(rz, pdv->z)) - gap;
    rgap = NA_REAL;
    if(pcost < 0.0) rgap = gap / (-pcost);
    if(dcost > 0.0) rgap = gap / dcost;
    pres = sqrt(resy * resy + resz * resz);
    dres = resx;
    if(i == 0){
      pres0 = std::max(1.0, pres);
      dres0 = std::max(1.0, dres);
    }
    pres = pres / pres0;
    dres = dres / dres0;
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
      ts = cList.smss(pdv->s).max();
      tz = cList.smss(pdv->z).max();
      state["pobj"] = pcost;
      state["dobj"] = dcost;
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
      cps->set_niter(i);
      cps->set_pdv(*pdv);
      cps->pdv.x.reshape(ne, 1); // removing variable 't'
      if((mnl == 1) && (cList.K == 1)){ // removing slack variables pertinent to 't'
	cps->pdv.s.set_size(0, 0);
	cps->pdv.z.set_size(0, 0);
	cps->set_sidx(umat());
      }
      if((mnl > 1) && (cList.K == 1)){ // removing slack variables pertinent to 't'
	cps->pdv.s.set_size(cList.dims[0] - 1, 1);
	cps->pdv.z.set_size(cList.dims[0] - 1, 1);
	cps->pdv.s = pdv->s.submat(1, 0, cList.dims[0] - 1, 0);
	cps->pdv.z = pdv->z.submat(1, 0, cList.dims[0] - 1, 0);
	umat sidxEpi = cList.sidx;
	sidxEpi.at(0, 1) -= 1;
	cps->set_sidx(sidxEpi);
      }
      if((mnl == 1) && (cList.K > 1)){ // removing slack variables pertinent to 't'
	cps->pdv.s.set_size(cList.G.n_rows - 1, 1);
	cps->pdv.z.set_size(cList.G.n_rows - 1, 1);
	cps->pdv.s = pdv->s.submat(1, 0, cList.G.n_rows - 1, 0);
	cps->pdv.z = pdv->z.submat(1, 0, cList.G.n_rows - 1, 0);
	umat sidxEpi = cList.sidx;
	sidxEpi.shed_row(0);
	sidxEpi -= 1;
	sidxEpi.at(0, 0) = 0;
	cps->set_sidx(sidxEpi);
      }
      if(trace){
	Rcpp::Rcout << "Optimal solution found." << std::endl;
      }
      return cps;
    }
    // Compute initial scalings
    if(i == 0){
      WList = cList.ntsc(pdv->s, pdv->z);
      Lambda = cList.getLambda(WList);
    }
    LambdaPrd = cList.sprd(Lambda, Lambda);
    LHS.submat(0, 0, ne - 1, ne - 1) = H;
    sigma = 0.0;
    // Finding solution of increments in two-round loop 
    // (same for affine and combined solution)
    for(int ii = 0; ii < 2; ii++){
      mu = gap / m;
      dpdv->s = -1.0 * LambdaPrd + OneE * sigma * mu;
      dpdv->x = -1.0 * rx;
      dpdv->y = -1.0 * ry;
      dpdv->z = -1.0 * rz;
      // Solving KKT-system
      try{
	dpdv->s = cList.sinv(dpdv->s, Lambda);
	Ws3 = cList.ssnt(dpdv->s, WList, false, true);
	dpdv->z = dpdv->z - Ws3;
	dpdv = sxyz(dpdv, LHS, WList); 
	dpdv->s = dpdv->s - dpdv->z;
      } catch(std::runtime_error &ex) {
	ts = cList.smss(pdv->s).max();
	tz = cList.smss(pdv->z).max();
	state["pobj"] = pcost;
	state["dobj"] = dcost;
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
      } catch(...) {
	::Rf_error("C++ exception (unknown reason)"); 
      }
      // Maximum step to boundary
      dpdv->s = cList.sslb(dpdv->s, Lambda, false);
      dpdv->z = cList.sslb(dpdv->z, Lambda, false); 
      ts = cList.smss(dpdv->s).max();
      tz = cList.smss(dpdv->z).max();
      ss << 0.0 << ts << tz << endr;
      tm = ss.max();
      if(tm == 0.0){
	step = 1.0;
      } else {
	step = std::min(1.0, sadj / tm);
      }
      // Backtracking until x is in the domain of f
      backTrack = true;
      while(backTrack){
	x = pdv->x + step * dpdv->x;
	for(int j = 0; j < mnl; j++){
	  Fval(j) = feval(x(span(0, ne - 1), span::all), nF[j]);
	}
	Fval[0] -= x.at(n - 1, 0); 
	if(is_finite(Fval)){
	  backTrack = false;
	} else {
	  step *= beta;
	}
      } // end while-loop domain of f
      if(ii == 0){
	sigma = pow((1.0 - step), 3.0);
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
  }  // end i-loop

  // Preparing result for non-convergence in maxiters iterations
  cps->set_pdv(*pdv);
  cps->pdv.x.reshape(ne, 1);    
  cps->set_sidx(cList.sidx);
  state["pobj"] = pobj(*pdv);
  state["dobj"] = dobj(*pdv);
  state["dgap"] = gap;
  state["certp"] = certp(*pdv);
  state["certd"] = certd(*pdv);
  ts = cList.smss(pdv->s).max();
  tz = cList.smss(pdv->z).max();
  state["pslack"] = -ts;
  state["dslack"] = -tz;
  if(!std::isnan(rgap)){
    state["rgap"] = rgap;
  }
  cps->set_state(state);
  cps->set_niter(maxiters);
  cps->set_status("unknown");
  if(trace){
    Rcpp::Rcout << "Optimal solution not determined in " << maxiters << " iteration(s)." << std::endl;
  }
  if((mnl == 1) && (cList.K == 1)){ // removing slack variables pertinent to 't'
    cps->pdv.s.set_size(0, 0);
    cps->pdv.z.set_size(0, 0);
    cps->set_sidx(umat());
  }
  if((mnl > 1) && (cList.K == 1)){ // removing slack variables pertinent to 't'
    cps->pdv.s.set_size(cList.dims[0] - 1, 1);
    cps->pdv.z.set_size(cList.dims[0] - 1, 1);
    cps->pdv.s = pdv->s.submat(1, 0, cList.dims[0] - 1, 0);
    cps->pdv.z = pdv->z.submat(1, 0, cList.dims[0] - 1, 0);
    umat sidxEpi = cList.sidx;
    sidxEpi.at(0, 1) -= 1;
    cps->set_sidx(sidxEpi);
  }
  if((mnl == 1) && (cList.K > 1)){ // removing slack variables pertinent to 't'
    cps->pdv.s.set_size(cList.G.n_rows - 1, 1);
    cps->pdv.z.set_size(cList.G.n_rows - 1, 1);
    cps->pdv.s = pdv->s.submat(1, 0, cList.G.n_rows - 1, 0);
    cps->pdv.z = pdv->z.submat(1, 0, cList.G.n_rows - 1, 0);
    umat sidxEpi = cList.sidx;
    sidxEpi.shed_row(0);
    sidxEpi -= 1;
    sidxEpi.at(0, 0) = 0;
    cps->set_sidx(sidxEpi);
  }

  return cps;
}
