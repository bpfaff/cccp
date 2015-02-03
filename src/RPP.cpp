/*
 * Function for solving risk parity portfolios (long-only)
*/
#include "cccp3.h"
using namespace arma;

CPS* rpp(mat x0, mat P, mat mrc, CTRL& ctrl){
  // Initializing objects
  int ne = P.n_cols;
  int n = ne + 1;
  // Constraints
  CONEC cList, cEpi;
  std::vector<std::string> cones; 
  cones.push_back("NLFC");
  cones.push_back("NNOC");
  cList.cone = cones;
  cList.G.eye(ne, ne);
  cList.G *= -1.0;
  cList.G.insert_cols(ne, 1);
  cList.G.insert_rows(0, 1);
  cList.G(0, ne) = -1.0;
  cList.h.zeros(n, 1);
  cList.dims << 1 << ne << endr;
  cList.sidx << 0 << 0 << endr
	     << 1 << ne << endr;
  cList.K = 2;
  cList.n = n;
  // Primal dual variables
  PDV* pdv = cList.initpdv(0);
  PDV* dpdv = cList.initpdv(0);
  pdv->x(span(0, ne - 1), span::all) = x0;
  // Solution
  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  // Objects used in iterations
  Rcpp::NumericVector state = cps->get_state();
  bool checkRgap = false, backTrack;
  int m = sum(cList.dims);
  double Fval, gap = m, resx, resz, pcost = 1.0, dcost = 1.0, rgap = NA_REAL, 
    pres = 1.0, dres = 1.0, pres0 = 1.0, dres0 = 1.0, sigma, mu, ts, tz, tm, step, a, x1;
  vec ss(3);
  mat H(ne, ne), LHS = zeros(ne, ne), rx, rz(cList.G.n_rows, 1), Lambda, LambdaPrd, Ws3, x, 
    ux(ne, 1), uz, RHS(ne, 1);
  mat OneE = cList.sone();
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
    // Setting f0 to first row of h-matrix
    cList.h(0, 0) = rpp_f0(pdv->x(span(0, ne - 1), span::all), P, mrc) - pdv->x.at(n - 1, 0);
    // Setting Df to first row of G-matrix
    cList.G(0, span(0, ne - 1)) = rpp_g0(pdv->x(span(0, ne - 1), span::all), P, mrc).t();
    // Computing Hessian
    H = pdv->z.at(0, 0) * rpp_h0(pdv->x(span(0, ne - 1), span::all), P, mrc);
    // Computing gap
    gap = sum(cList.sdot(pdv->s, pdv->z));
    // Computing residuals
    // Dual Residuals
    rx = cList.G.t() * pdv->z;
    rx.at(rx.n_rows - 1, 0) += 1.0;
    resx = norm(rx);
    // Central Residuals 
    rz.at(0, 0) = pdv->s.at(0, 0) + cList.h.at(0, 0);
    rz(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) = 
	pdv->s(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) +
	cList.G(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all) * pdv->x -
	cList.h(span(cList.sidx.at(1, 0), cList.sidx.at(1, 1)), span::all);
    resz = cList.snrm2(rz);
    // Statistics for stopping criteria
    pcost = pdv->x.at(pdv->x.n_rows - 1, 0);
    dcost = pcost + sum(cList.sdot(rz, pdv->z)) - gap;
    rgap = NA_REAL;
    if(pcost < 0.0) rgap = gap / (-pcost);
    if(dcost > 0.0) rgap = gap / dcost;
    pres = sqrt(resz * resz);
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
      cps->pdv.x = cps->pdv.x / accu(cps->pdv.x); // Budget constraint
      cps->pdv.s.set_size(cList.G.n_rows - 1, 1);
      cps->pdv.z.set_size(cList.G.n_rows - 1, 1);
      cps->pdv.s = pdv->s.submat(1, 0, cList.G.n_rows - 1, 0);
      cps->pdv.z = pdv->z.submat(1, 0, cList.G.n_rows - 1, 0);
      umat sidxEpi = cList.sidx;
      sidxEpi.shed_row(0);
      sidxEpi -= 1;
      sidxEpi.at(0, 0) = 0;
      cps->set_sidx(sidxEpi);
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
    // Solution step 1 in two-round loop 
    // (same for affine and combined solution)
    for(int ii = 0; ii < 2; ii++){
      mu = gap / m;
      dpdv->s = -1.0 * LambdaPrd + OneE * sigma * mu;
      dpdv->x = -1.0 * rx;
      dpdv->z = -1.0 * rz;
      // Solving KKT-system
      try{
	dpdv->s = cList.sinv(dpdv->s, Lambda);
	Ws3 = cList.ssnt(dpdv->s, WList, false, true);
	dpdv->z = dpdv->z - Ws3;
	// Solving reduced reduced system
	cEpi = cList;
	WEpi = WList; 
	a = dpdv->z.at(0, 0); // Slack with respect to f0
	x1 = dpdv->x.at(dpdv->x.n_rows, 0); // Epigraph-variable 't'
	uz = dpdv->z;
	cEpi.n -= 1;
	cEpi.G.set_size(cList.G.n_rows, ne);
	cEpi.G = cList.G(span::all, span(0, ne - 1)); // removing last column
	ux = dpdv->x(span(0, ne - 1), span::all);
	ux = ux + dpdv->x.at(n - 1, 0) * cEpi.G.row(0).t();
	WEpi.erase(WEpi.begin());
	cEpi.K -= 1;
	cEpi.G = cEpi.G(span(1, cEpi.G.n_rows - 1), span::all); // removing first row pertinent to f0
	cEpi.sidx = cEpi.sidx(span(1, cEpi.sidx.n_rows - 1), span::all);
	cEpi.sidx -= 1;
	cEpi.sidx.at(0, 0) = 0;
	cEpi.cone.erase(cEpi.cone.begin());
	cEpi.dims.shed_row(0);
	LHS = H + cEpi.gwwg(WEpi);
	uz = uz(span(1, uz.n_rows - 1), span::all); 
	RHS = ux + cEpi.gwwz(WEpi, uz);
	dpdv->x.submat(0, 0, ne - 1, 0) = solve(LHS, RHS);
	// Preparing dpdv
	uz = cEpi.G * dpdv->x.submat(0, 0, ne - 1, 0) - uz;
	dpdv->z(span(1, dpdv->z.n_rows - 1), span::all) = cEpi.ssnt(uz, WEpi, true, true);
	dpdv->z.at(0, 0) = -dpdv->x.at(dpdv->x.n_rows - 1, 0) * WList[0]["dnl"].at(0, 0);
	x1 = dot(cEpi.G(0, span::all), dpdv->x.submat(0, 0, ne - 1, 0)) + 
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
	Fval = rpp_f0(x(span(0, ne - 1), span::all), P, mrc) - x.at(n - 1, 0);
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
  cps->pdv.x = cps->pdv.x / accu(cps->pdv.x); // Budget constraint
  cps->pdv.s.set_size(cList.G.n_rows - 1, 1);
  cps->pdv.z.set_size(cList.G.n_rows - 1, 1);
  cps->pdv.s = pdv->s.submat(1, 0, cList.G.n_rows - 1, 0);
  cps->pdv.z = pdv->z.submat(1, 0, cList.G.n_rows - 1, 0);
  umat sidxEpi = cList.sidx;
  sidxEpi.shed_row(0);
  sidxEpi -= 1;
  sidxEpi.at(0, 0) = 0;
  cps->set_sidx(sidxEpi);

  return cps;
}
