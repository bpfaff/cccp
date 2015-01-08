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
