/* 
 *
 * Header file for package cccp3
 *
*/
#ifndef CCCP3_H
#define CCCP3_H

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

double sdot_nls(arma::mat s, arma::mat z);
double sdot_p(arma::mat s, arma::mat z, int m);
double snrm2_nls(arma::mat s);
double snrm2_p(arma::mat s, int m);

#include "NLF.h"
#include "NNO.h"
#include "SOC.h"
#include "PSD.h"
#include "CPG.h"

#endif
