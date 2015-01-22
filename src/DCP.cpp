#include "cccp3.h"
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
  double ans = pdv.x(pdv.x.n_rows, 0);
  return ans;
}
/*
Dual objective
*/
double DCP::dobj(PDV& pdv){
  double term1 = 0.0, term2 = 0.0, term3 = 0.0, ans;
  term1 = pdv.x(pdv.x.n_rows, 0);
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
  mat ans(p,1);
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

  if(cList.K > 0){
    Gz = cList.G.t() * pdv.z;
  }
  if(A.n_rows > 0){
    Ay = A.t() * pdv.y;
  }
  ans = Gz + Ay;
  ans(ans.n_rows, 0) += 1.0;

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
  int m = cList.G.n_rows;
  int mnl = cList.dims(0);
  int K = cList.K;
  double a = pdv->z(0, 0); // Slack with respect to f0
  double x1 = pdv->x(pdv->x.n_rows, 0); // Epigraph-variable 't'
  mat ux = pdv->x, uz = pdv->z, RHS(ne + A.n_rows, 1), ans(ne + A.n_rows, 1);
  CONEC cEpi = cList;
  std::vector<std::map<std::string,mat> > WEpi = WList;
  
  // Distinguishing four cases:
  // mnl == 1 and K == 1: only f0 and no cone constraints
  // mnl == 1 and K > 1 : only f0 and cone constraints
  // mnl > 1 and K == 1 : f0 and other nonlinear constraints; no cone constraints
  // mnl > 1 and K > 1 : f0 and other nonlinear constraints and cone constraints
  // Problem to be solved is reduced to x0

  cEpi.n -= 1;
  cEpi.G = cEpi.G(span::all, span(0, cEpi.G.n_cols - 1)); // removing last column
  // mnl == 1 and K == 1
  if((mnl == 1) && (K == 1)){
    // upper LHS is only Hessian
    ux = ux(span(0, ux.n_rows - 1), span::all);
    ux += pdv->x(pdv->x.n_rows, 0) * cEpi.G.row(0).t();
    RHS.submat(0, 0, ne - 1, 0) = ux;
  }
  // mnl == 1 and K > 1
  if((mnl == 1) && (K > 1)){
    WEpi.erase(WEpi.begin());
    cEpi.K -= 1;
    cEpi.G = cEpi.G(span(1, cEpi.G.n_rows), span::all); // removing first row pertinent to f0
    cEpi.sidx = cEpi.sidx(span(1, cEpi.sidx.n_rows), span::all);
    cEpi.cone.erase(cEpi.cone.begin());
    cEpi.dims = cEpi.dims(1, cEpi.dims.size());
    LHS.submat(0, 0, ne - 1, ne - 1) += cEpi.gwwg(WEpi);
    uz = uz(span(1, uz.n_rows), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  // mnl > 1 and K == 1
  if((mnl > 1) && (K == 1)){
    WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
    WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
    cEpi.dims(0) -= 1;
    cEpi.sidx(0, 1) -= 1;
    LHS.submat(0, 0, ne - 1, ne - 1) += cEpi.gwwg(WEpi);
    uz = uz(span(1, uz.n_rows), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  // mnl > 1 and K > 1
  if((mnl > 1) && (K >1)){
    WEpi[0]["dnl"] = WEpi[0]["dnl"](span(1, mnl - 1), span::all);
    WEpi[0]["dnli"] = WEpi[0]["dnli"](span(1, mnl - 1), span::all);
    cEpi.dims(0) -= 1;
    cEpi.sidx = cEpi.sidx - 1;
    cEpi.sidx(0, 0) = 0;
    uz = uz(span(1, uz.n_rows), span::all);
    RHS.submat(0, 0, ne - 1, 0) = ux + cEpi.gwwz(WEpi, uz);
  }
  if(pdv->y.n_rows > 0){
    RHS.submat(n, 0, RHS.n_rows - 1, 0) = pdv->y;
  }
  // Solving KKT-system
  ans = solve(LHS, RHS);
  // Preparing pdv
  pdv->x.submat(0, 0, ne - 1, 0) = ans.submat(0, 0, ne - 1, 0);
  if(pdv->y.n_rows > 0){
    pdv->y = ans.submat(ne, 0, RHS.n_rows - 1, 0);
  }
  uz = cEpi.G * pdv->x.submat(0, 0, ne - 1, 0) - uz;
  pdv->z(span(1, pdv->z.n_rows), span::all) = cEpi.ssnt(uz, WEpi, true, true);
  pdv->z(0, 0) = -pdv->x(pdv->x.n_rows, 0) * WList[0]["dnl"](0, 0);
  x1 = dot(cList.G(0, span(0, ne - 1)), pdv->x.submat(0, 0, ne - 1, 0)) + 
    pow(WList[0]["dnl"](0, 0), 2) + x1 - a; 
  pdv->x(pdv->x.n_rows, 0) = x1;

  return pdv;
}
