/* 
 *
 * Header file for package cccp3
 *
*/

#ifndef CCCP3_H
#define CCCP3_H

#include <RcppCommon.h>
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(NLFV)
RCPP_EXPOSED_CLASS(NLFS)
RCPP_EXPOSED_CLASS(NNOV)
RCPP_EXPOSED_CLASS(NNOS)
RCPP_EXPOSED_CLASS(SOCV)
RCPP_EXPOSED_CLASS(SOCS)
RCPP_EXPOSED_CLASS(PSDV)

#include <RcppArmadillo.h>

double udot_n(NLFV* s, NLFV* z); // for Rcpp Module 
double udot(NLFV s);
double udot(NLFV s, NLFV z);

double udot_l(NNOV* s, NNOV* z); // for Rcpp Module 
double udot(NNOV s);
double udot(NNOV s, NNOV z);

double udot_s(SOCV* s, SOCV* z); // for Rcpp Module 
double udot(SOCV s);
double udot(SOCV s, SOCV z);

double udot_p(PSDV* s, PSDV* z); // for Rcpp Module 
double udot(PSDV s);
double udot(PSDV s, PSDV z);

double jdot_s(SOCV* s, SOCV* z); // for Rcpp Module 
double jdot(SOCV s);
double jdot(SOCV s, SOCV z);

NLFV uone_n(NLFV* s); // for Rcpp Module
NNOV uone_l(NNOV* s); // for Rcpp Module
SOCV uone_s(SOCV* s); // for Rcpp Module
PSDV uone_p(PSDV* s); // for Rcpp Module

NLFV uprd_n(NLFV* s, NLFV* z); // for Rcpp Module
NNOV uprd_l(NNOV* s, NNOV* z); // for Rcpp Module
SOCV uprd_s(SOCV* s, SOCV* z); // for Rcpp Module
PSDV uprd_p(PSDV* s, PSDV* z); // for Rcpp Module

NLFV uinv_n(NLFV* s, NLFV* z); // for Rcpp Module
NNOV uinv_l(NNOV* s, NNOV* z); // for Rcpp Module
SOCV uinv_s(SOCV* s, SOCV* z); // for Rcpp Module
PSDV uinv_p(PSDV* s, PSDV* z); // for Rcpp Module

SEXP umss_n(NLFV* s); // for Rcpp Module
SEXP umss_l(NNOV* s); // for Rcpp Module
SEXP umss_s(SOCV* s); // for Rcpp Module
SEXP umss_p(PSDV* s); // for Rcpp Module

NLFV umsa_n(NLFV* s, double alpha, bool init); // for Rcpp Module
NNOV umsa_l(NNOV* s, double alpha, bool init); // for Rcpp Module
SOCV umsa_s(SOCV* s, double alpha, bool init); // for Rcpp Module
PSDV umsa_p1(PSDV* s, double alpha); // for Rcpp Module 
PSDV umsa_p2(PSDV* s, double alpha, arma::vec sigma, PSDV* lmbda); // for Rcpp Module 

NLFS ntsc_n(NLFV* s, NLFV* z); // for Rcpp Module 
NNOS ntsc_l(NNOV* s, NNOV* z); // for Rcpp Module 
SOCS ntsc_s(SOCV* s, SOCV* z); // for Rcpp Module 

#include "NLF.h"
#include "NNO.h"
#include "SOC.h"
#include "PSD.h"

#endif
