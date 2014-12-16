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

double sdot_nlp(arma::mat s, arma::mat z);
double sdot_s(arma::mat s, arma::mat z, int m);
double jdot_p(arma::mat s, arma::mat z);
double snrm2_nlp(arma::mat s);
double snrm2_s(arma::mat s, int m);
double jnrm2_p(arma::mat s);
arma::mat sprd_nl(arma::mat s, arma::mat z);
arma::mat sprd_p(arma::mat s, arma::mat z);
arma::mat sprd_s(arma::mat s, arma::mat z, int m);
arma::mat sone_nl(arma::mat s);
arma::mat sone_p(arma::mat s);
arma::mat sone_s(int m);
arma::mat sinv_nl(arma::mat s, arma::mat z);
arma::mat sinv_p(arma::mat s, arma::mat z);
arma::mat sinv_s(arma::mat s, arma::mat z, int m);
double smss_nl(arma::mat s);
double smss_p(arma::mat s);
double smss_s(arma::mat s, int m);
arma::mat sams1_nl(arma::mat s, double alpha);
arma::mat sams1_p(arma::mat s, double alpha);
arma::mat sams1_s(arma::mat s, double alpha, int m);
arma::mat sams2_nl(arma::mat s, double alpha);
arma::mat sams2_p(arma::mat s, double alpha);
arma::mat sams2_s(arma::mat s, double alpha, arma::mat lambda, arma::vec sigma, int m);

#include "NLF.h"
#include "NNO.h"
#include "SOC.h"
#include "PSD.h"
#include "CPG.h"

#endif
