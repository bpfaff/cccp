#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(PDV)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition for primal/dual variables
*/
class PDV {
 public:

  // constructors
 PDV() : x(arma::vec()), y(arma::vec()), s(Rcpp::List::create()), z(Rcpp::List::create()), kappa(1.0), tau(1.0) {}
 PDV(arma::vec x_, arma::vec y_, Rcpp::List s_, Rcpp::List z_, double kappa_, double tau_):  x(x_), y(y_), s(s_), z(z_), kappa(kappa_), tau(tau_) {}
  // members
  arma::vec get_x() {return x;}
  void set_x(arma::vec x_) {x = x_;}
  arma::vec get_y() {return y;}
  void set_y(arma::vec y_) {y = y_;}
  Rcpp::List get_s() {return s;}
  void set_s(Rcpp::List s_) {s = s_;}
  Rcpp::List get_z() {return z;}
  void set_z(Rcpp::List z_) {z = z_;}
  double get_kappa() {return kappa;}
  void set_kappa(double kappa_) {kappa = kappa_;}
  double get_tau() {return tau;}
  void set_tau(double tau_) {tau = tau_;}

 private:
  arma::vec x;
  arma::vec y;
  Rcpp::List s;
  Rcpp::List z;
  double kappa;
  double tau;
};

