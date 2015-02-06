#include "cccp.h"
/*
 *
 * Methods for Linear Programs with nonlinear constraints
 *
*/
using namespace arma;
/*
Primal objective
*/
double DNL::pobj(PDV& pdv){
  double ans = dot(pdv.x, q);
  return ans;
}
/*
Dual objective
*/
double DNL::dobj(PDV& pdv){
  double term1 = 0.0, term2 = 0.0, term3 = 0.0, ans;
  term1 = dot(pdv.x, q);
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
mat DNL::rprim(PDV& pdv){
  int p = A.n_rows;
  mat ans(p,1);
  ans.zeros();

  ans = b - A * pdv.x;

  return ans;
}
/*
Centrality Residuals
*/
mat DNL::rcent(PDV& pdv){
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
mat DNL::rdual(PDV& pdv){
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
double DNL::certp(PDV& pdv){
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
double DNL::certd(PDV& pdv){
  double nomin, denom, ans;

  nomin = norm(rdual(pdv));
  denom = std::max(1.0, norm(q));
  ans = nomin / denom;

  return ans;
}
/*
  Solving 'KKT-System'
*/
PDV* DNL::sxyz(PDV* pdv, mat LHS, mat RHS, std::vector<std::map<std::string,mat> > WList){
  int n = q.n_rows;
  mat rhs1, ans;

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
  Main routine for solving a Linear Program with nonlinear constraints
*/
CPS* DNL::cps(CTRL& ctrl){
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
  int m = sum(cList.dims), mnl = cList.dims(0), n = cList.n, sizeLHS = A.n_rows + A.n_cols;
  double gap = m, resx, resy, resz, pcost, dcost, rgap = NA_REAL, 
    pres, dres, pres0 = 1.0, dres0 = 1.0, sigma, mu, ts, tz, tm, step;
  vec ss(3), Fval(mnl);
  mat H = zeros(n, n), rx, ry, rz, Lambda, LambdaPrd, Ws3, x;
  mat OneE = cList.sone();
  mat LHS(sizeLHS, sizeLHS);
  // Initialising LHS matrices
  LHS.zeros();
  if(A.n_rows > 0){ // equality constraints
    LHS.submat(n, 0, sizeLHS-1, n-1) = A;
    LHS.submat(0, n, n-1, sizeLHS-1) = A.t();
  }
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
  double beta = Rcpp::as<double>(params["beta"]);
  //
  // Starting iterations
  //
  for(int i = 0; i < maxiters; i++){
    H.zeros();
    for(int j = 0; j < mnl; j++){
      // Setting f to first mnl-rows of h-matrix
      cList.h(j, 0) = feval(pdv->x, nF[j]);
      // Setting Df to first mnl-rows of G-matrix
      cList.G.row(j) = geval(pdv->x, gF[j]).t();
      // Computing Hessian
      H += pdv->z.at(j, 0) * heval(pdv->x, hF[j]);
    }
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
    LHS.submat(0, 0, n-1, n-1) = H + cList.gwwg(WList);
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
	dpdv = sxyz(dpdv, LHS, RHS, WList); 
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
	  Fval(j) = feval(x, nF[j]);
	}
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

  return cps;
}
