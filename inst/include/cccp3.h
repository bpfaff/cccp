/* 
 *
 * Header file for package cccp3
 *
*/

#ifndef CCCP3_H
#define CCCP3_H

#include <RcppCommon.h>
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(PSDV)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

double udot_p(PSDV* s, PSDV* z); // for Rcpp Module 
PSDV uone_p(PSDV* s); // for Rcpp Module
PSDV uprd_p(PSDV* s, PSDV* z); // for Rcpp Module
PSDV uinv_p(PSDV* s, PSDV* z); // for Rcpp Module
SEXP umss_p(PSDV* s); // for Rcpp Module
PSDV umsa_p1(PSDV* s, double alpha); // for Rcpp Module 
PSDV umsa_p2(PSDV* s, double alpha, arma::vec sigma, PSDV* lmbda); // for Rcpp Module 

#include "NLF.h"
#include "NNO.h"
#include "SOC.h"
#include "PSD.h"

#endif
