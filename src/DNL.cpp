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
double DNL::pobj(PDV& pdv){
  double ans = dot(pdv.x, q);
  return ans;
}
/*
Dual objective
*/
double DNL::dobj(PDV& pdv){
  double term1 = 0.0, term2 = 0.0, ans;
  term1 = dot(b, pdv.y);
  term2 = sum(cList.sdot(pdv.z, cList.h));
  ans = -term1 - term2;
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

  ans = pdv.s + cList.G * pdv.x - cList.h;

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
