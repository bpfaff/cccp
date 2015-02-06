/*
 * Function for solving a geometric program
*/
#include "cccp.h"
using namespace arma;

CPS* gpp(std::vector<mat> FList, std::vector<mat> gList, CONEC& cList, mat A, mat b, CTRL& ctrl){
  // Initializing objects
  int n = cList.n, ne = n - 1, m = sum(cList.dims), 
    mnl = cList.dims(0), sizeLHS = A.n_rows + A.n_cols - 1;
  // Constraints
  CONEC cEpi;
  // Primal dual variables
  PDV* pdv = cList.initpdv(A.n_rows);
  PDV* dpdv = cList.initpdv(A.n_rows);
  // Solution
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  // Objects used in iterations
  Rcpp::NumericVector state = cps->get_state();
  bool checkRgap = false, backTrack;
  double gap = m, resx, resy, resz, pcost = 1.0, dcost = 1.0, rgap = NA_REAL, 
    pres = 1.0, dres = 1.0, pres0 = 1.0, dres0 = 1.0, sigma, mu, ts, tz, tm, step, a, x1,
    ymax, ysum;
  vec ss(3);
  mat H(ne, ne), rx, ry, rz(cList.G.n_rows, 1), Lambda, LambdaPrd, Ws3, x, ans(sizeLHS, 1),
    ux(ne, 1), uz, RHS(sizeLHS, 1), y, Fval(mnl, 1);
  mat OneE = cList.sone();
  // Initialising LHS matrices
  mat LHS = zeros(sizeLHS, sizeLHS);
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(ne, 0, sizeLHS - 1, ne - 1) = A(span::all, span(0, ne - 1));
    LHS.submat(0, ne, ne - 1, sizeLHS - 1) = A(span::all, span(0, ne - 1)).t();
  }
  std::vector<mat> FGP;
  std::vector<std::map<std::string,mat> > WList, WEpi;
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
      FGP = fgp(pdv->x(span(0, ne - 1), span::all), FList[j], gList[j]);
      // Setting f to first mnl-rows of h-matrix
      cList.h(j, 0) = FGP[0].at(0, 0);
      // Setting Df to first mnl-rows of G-matrix
      cList.G(j, span(0, ne - 1)) = FGP[1].t();
      // Computing Hessian
      H += pdv->z.at(j, 0) * FGP[2];
    }
    cList.h(0, 0) = cList.h(0, 0) - pdv->x.at(n - 1, 0);
    // Computing gap
    gap = sum(cList.sdot(pdv->s, pdv->z));
    // Computing residuals
    // Dual Residuals
    rx = cList.G.t() * pdv->z + A.t() * pdv->y;
    rx.at(rx.n_rows - 1, 0) += 1.0; // last element of 'q' is 1.0 * t
    resx = norm(rx);
    // Primal Residuals
    ry = b - A * pdv->x;
    resy = norm(ry);
    // Central Residuals 
    rz(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all) =
      pdv->s(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all) + 
      cList.h(span(cList.sidx.at(0, 0), cList.sidx.at(0, 1)), span::all);
    if(cList.K > 1){
      rz(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) =
	pdv->s(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) + 
	cList.G(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) * pdv->x -
	cList.h(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all);
    }
    resz = cList.snrm2(rz);
    // Statistics for stopping criteria
    pcost = pdv->x.at(n - 1, 0);
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
      cps->pdv.x = exp(cps->pdv.x);
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
      cps->pdv.s = exp(cps->pdv.s);
      cps->pdv.z = exp(cps->pdv.z);
      if(A.n_rows > 0){
	cps->pdv.y = exp(cps->pdv.y);
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
    sigma = 0.0;
    //
    // Creating objects for epigraph form
    //
    cEpi = cList;
    WEpi = WList; 
    cEpi.n -= 1;
    cEpi.G.shed_col(n - 1);
    cEpi.G.shed_row(0);
    cEpi.h.shed_row(0);
    // Distinguishing three cases:
    // mnl == 1 and K > 1 : only f0 and NNO constraints
    // mnl > 1 and K == 1 : f0 and posinomial constraints; no NNO constraints
    // mnl > 1 and K > 1 : f0, posinomial constraints and cone constraints
    // Problem to be solved is reduced to x0

    // mnl == 1 and K > 1
    if((mnl == 1) && (cList.K > 1)){
      WEpi.erase(WEpi.begin());
      cEpi.K -= 1;
      cEpi.sidx.shed_row(0);
      cEpi.sidx -= 1;
      cEpi.sidx.at(0, 0) = 0;
      cEpi.cone.erase(cEpi.cone.begin());
      cEpi.dims.shed_row(0);
    }
    // mnl > 1 and K == 1
    if((mnl > 1) && (cList.K == 1)){
      WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
      WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
      cEpi.dims(0) -= 1;
      cEpi.sidx(0, 1) -= 1;
    }
    // mnl > 1 and K > 1
    if((mnl > 1) && (cList.K > 1)){
      WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
      WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
      cEpi.dims(0) -= 1;
      cEpi.sidx = cEpi.sidx - 1;
      cEpi.sidx(0, 0) = 0;
    }
    LHS.submat(0, 0, ne - 1, ne - 1) = H + cEpi.gwwg(WEpi);	
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
	// Solving reduced system
	a = dpdv->z.at(0, 0); // Slack with respect to f0
	x1 = dpdv->x.at(n - 1, 0); // Epigraph-variable 't'
	ux = dpdv->x(span(0, ne - 1), span::all);
	ux = ux + dpdv->x.at(n - 1, 0) * cList.G.submat(0, 0, 0, ne - 1).t();
	uz = dpdv->z(span(1, dpdv->z.n_rows - 1), span::all);
	RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
	if(pdv->y.n_rows > 0){
	  RHS.submat(ne, 0, RHS.n_rows - 1, 0) = dpdv->y;
	}
	// Solving KKT-system
	ans = solve(LHS, RHS);
	dpdv->x.submat(0, 0, ne - 1, 0) = ans.submat(0, 0, ne - 1, 0);
	if(dpdv->y.n_rows > 0){
	  dpdv->y = ans.submat(ne, 0, RHS.n_rows - 1, 0);
	}
	// Preparing dpdv
	uz = cEpi.G * dpdv->x.submat(0, 0, ne - 1, 0) - uz;
	dpdv->z(span(1, dpdv->z.n_rows - 1), span::all) = cEpi.ssnt(uz, WEpi, true, true);
	dpdv->z.at(0, 0) = -dpdv->x.at(dpdv->x.n_rows - 1, 0) * WList[0]["dnl"].at(0, 0);
	x1 = dot(cList.G.submat(0, 0, 0, ne - 1), dpdv->x.submat(0, 0, ne - 1, 0)) + 
	  pow(WList[0]["dnl"].at(0, 0), 2) * dpdv->x.at(n - 1, 0) - a;
	dpdv->x.at(n - 1, 0) = x1;
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
	  y = FList[j] * x(span(0, ne - 1), span::all) + gList[j];
	  ymax = y.max();
	  y = exp(y - ymax);
	  ysum = norm(y, 1);
	  Fval.at(j, 0) = ymax + log(ysum);
	}
	Fval.at(0, 0) -= x.at(n - 1, 0); 
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
  } // end i-loop

  // Preparing result for non-convergence in maxiters iterations
  cps->set_pdv(*pdv);
  cps->pdv.x.reshape(ne, 1);    
  cps->pdv.x = exp(cps->pdv.x);
  cps->set_sidx(cList.sidx);
  state["pobj"] = pcost;
  state["dobj"] = dcost;
  state["dgap"] = gap;
  state["certp"] = pres;
  state["certd"] = dres;
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
  cps->pdv.s.set_size(cList.G.n_rows - 1, 1);
  cps->pdv.z.set_size(cList.G.n_rows - 1, 1);
  cps->pdv.s = pdv->s.submat(1, 0, cList.G.n_rows - 1, 0);
  cps->pdv.z = pdv->z.submat(1, 0, cList.G.n_rows - 1, 0);
  cps->pdv.s = exp(cps->pdv.s);
  cps->pdv.z = exp(cps->pdv.z);
  if(A.n_rows > 0){
    cps->pdv.y = exp(cps->pdv.y);
  }
  umat sidxEpi = cList.sidx;
  sidxEpi -= 1;
  sidxEpi.at(0, 0) = 0;
  cps->set_sidx(sidxEpi);

  return cps;
}
