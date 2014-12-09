#include "cccp3.h"
/*
 * Methods / functions for cone variables
*/

/*
 * umsa: Applying maximum step-size to vectors in S
// PSD-variables (for initial max step length)
PSDV umsa_p1(PSDV* s, double alpha){ 
  s->u.reshape(s->dims, s->dims);
  s->u.diag() = s->u.diag() + (1 + alpha);
  s->u.reshape(s->dims * s->dims, 1);

  return *s;
}
*/

/*
// PSD-variables (for max step length)
PSDV umsa_p2(PSDV* s, double alpha, arma::vec sigma, PSDV* lmbda){
  s->u.reshape(s->dims, s->dims);
  lmbda->u.reshape(lmbda->dims, lmbda->dims);
  for(int i = 0; i < s->dims; i++){
    sigma(i) = 1 + alpha * sigma(i);
    s->u.col(i) = s->u.col(i) * sqrt(sigma(i) / lmbda->u(i, i)); 
  }
  s->u.reshape(s->dims * s->dims, 1);
  
  return *s;
}
*/
