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
Initializing PDV
*/
PDV* DLP::initpdv(){
  PDV* pdv = new PDV();
  mat s(cList.G.n_rows, 1);
  mat ans;
  int n = q.n_cols;

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
