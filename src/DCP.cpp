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
