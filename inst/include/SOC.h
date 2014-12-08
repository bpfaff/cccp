#include <RcppCommon.h>
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(SOCV)
RCPP_EXPOSED_CLASS(SOCS)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition and methods for second-order cone constraints
*/
// Class definition for vectors/variables
class SOCV {
public:

  // constructors
  SOCV() : u(arma::mat()), dims(0L) {}
  SOCV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
  SOCV(arma::mat u_, int dims_): u(u_), dims(dims_) {}
  // members
  arma::mat get_u() {return u;}
  void set_u(arma::mat u_) {u = u_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

  SOCV* uone() const;
  double umss() const;
  double udot(const SOCV& z) const;
  double jdot(const SOCV& z) const;
  SOCV* uprd(const SOCV& z) const;
  SOCV* uinv(const SOCV& z) const;
  SOCV* umsa(double alpha, bool init) const;
  SOCS* ntsc(const SOCV& z) const;

 private:
  arma::mat u;	
  int dims;
};

// Class definition for NT-scaling and Lagrange-Multipliers
class SOCS {
public:

  // constructors
  SOCS(arma::mat v_, double beta_, SOCV lambda_): v(v_), beta(beta_), lambda(lambda_) {}

  arma::mat get_v() {return v;}
  void set_v(arma::mat v_) {v = v_;}
  double get_beta() {return beta;}
  void set_beta(double beta_) {beta = beta_;}
  SOCV get_lambda() {return lambda;}
  void set_lambda(SOCV lambda_) {lambda = lambda_;}

 private:
  arma::mat v;
  double beta;
  SOCV lambda;
};
