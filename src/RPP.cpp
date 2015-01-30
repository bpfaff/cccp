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
  CONEC cList;
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

  CPS* cps = new CPS();
  cps->set_pdv(*pdv);
  cps->set_sidx(cList.sidx);
  Rcpp::NumericVector state = cps->get_state();
  bool checkRgap = false, backTrack;
  int m = sum(cList.dims), mnl = 1, sizeLHS = ne;
  double Fval, gap = m, resx, resy, resz, pcost, dcost, rgap = NA_REAL, 
    pres, dres, pres0 = 1.0, dres0 = 1.0, sigma, mu, ts, tz, tm, step;
  vec ss(3);
  mat LHS = zeros(ne, ne), rx, ry, rz, Lambda, LambdaPrd, Ws3, x;
  mat OneE = cList.sone();
  // Initialising LHS matrices
  LHS.zeros();
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
    // Setting f0 to first row of h-matrix
    cList.h(0, 0) = rpp_f0(pdv->x(span(0, ne - 1), span::all), P, mrc);
    cList.h(0, 0) -= pdv->x.at(n - 1, 0);
    // Setting Df to first row of G-matrix
    cList.G(0, span(0, ne - 1)) = rpp_g0(pdv->x(span(0, ne - 1), span::all), P, mrc).t();
    // Computing Hessian
    LHS = pdv->z.at(0, 0) * rpp_h0(pdv->x(span(0, ne - 1), span::all), P, mrc);
    // Computing gap
    gap = sum(cList.sdot(pdv->s, pdv->z));
  }

  return cps;
}
