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
using namespace arma;

double sdot_nlp(mat s, mat z);
double sdot_s(mat s, mat z, int m);
double jdot_p(mat s, mat z);
double snrm2_nlp(mat s);
double snrm2_s(mat s, int m);
double jnrm2_p(mat s);
mat sprd_nl(mat s, mat z);
mat sprd_p(mat s, mat z);
mat sprd_s(mat s, mat z, int m);
mat sone_nl(int m);
mat sone_p(int m);
mat sone_s(int m);
mat sinv_nl(mat s, mat z);
mat sinv_p(mat s, mat z);
mat sinv_s(mat s, mat z, int m);
double smss_nl(mat s);
double smss_p(mat s);
double smss_s(mat s, int m);
mat sams1_nl(mat s, double alpha);
mat sams1_p(mat s, double alpha);
mat sams1_s(mat s, double alpha, int m);
mat sams2_nl(mat s, double alpha);
mat sams2_p(mat s, double alpha);
mat sams2_s(mat s, double alpha, mat lambda, vec sigma, int m);
std::map<std::string,mat> ntsc_n(mat s, mat z);
std::map<std::string,mat> ntsc_l(mat s, mat z);
std::map<std::string,mat> ntsc_p(mat s, mat z);
std::map<std::string,mat> ntsc_s(mat s, mat z, int m);
std::map<std::string,mat> ntsu_n(std::map<std::string,mat> W, mat s, mat z);
std::map<std::string,mat> ntsu_l(std::map<std::string,mat> W, mat s, mat z);
std::map<std::string,mat> ntsu_p(std::map<std::string,mat> W, mat s, mat z);
std::map<std::string,mat> ntsu_s(std::map<std::string,mat> W, mat s, mat z, int m);
mat sslb_nl(mat s, mat lambda, bool invers);
mat sslb_p(mat s, mat lambda, bool invers);
mat sslb_s(mat s, mat lambda, bool invers, int m);
mat ssnt_n(mat s, std::map<std::string,mat> W, bool invers);
mat ssnt_l(mat s, std::map<std::string,mat> W, bool invers);
mat ssnt_p(mat s, std::map<std::string,mat> W, bool invers);
mat ssnt_s(mat s, std::map<std::string,mat> W, bool invers, bool transp);
double feval(mat x, Rcpp::Function Rf);
vec geval(mat x, Rcpp::Function Rf);
mat heval(mat x, Rcpp::Function Rf);

// Objective, Gradient and Hessian functions for risk parity
double rpp_f0(mat x, mat P, mat mrc);
mat rpp_g0(mat x, mat P, mat mrc);
mat rpp_h0(mat x, mat P, mat mrc);

// Function value, Gradient and Hessian for geometric programs
std::vector<mat> fgp(mat x, mat F, mat g);

#include "CPG.h"

#endif
