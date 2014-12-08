#include <RcppCommon.h>
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(NLFV)
RCPP_EXPOSED_CLASS(NLFS)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition and methods for nonlinear constraints
*/
// Class definition for vectors/variables 
class NLFV {
 public:

  // constructors
 NLFV() : u(arma::mat()), dims(0L) {}
 NLFV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
 NLFV(arma::mat u_, int dims_): u(u_), dims(dims_) {}
  // members
  arma::mat get_u() {return u;}
  void set_u(arma::mat u_) {u = u_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

  NLFV* uone() const;
  double umss() const;
  double udot(const NLFV& z) const;
  NLFV* uprd(const NLFV& z) const;
  NLFV* uinv(const NLFV& z) const;
  NLFV* umsa(double alpha, bool init) const;
  NLFS* ntsc(const NLFV& z) const;

 private:
  arma::mat u;
  int dims;
};


// Class definition for NT-scaling and Lagrange-Multipliers
class NLFS {
 public:

  // constructors
 NLFS(arma::mat dnl_, arma::mat dnli_, NLFV lambda_): dnl(dnl_), dnli(dnli_), lambda(lambda_) {}

  arma::mat get_dnl() {return dnl;}
  void set_dnl(arma::mat dnl_) {dnl = dnl_;}
  arma::mat get_dnli() {return dnli;}
  void set_dnli(arma::mat dnli_) {dnli = dnli_;}
  NLFV get_lambda() {return lambda;}
  void set_lambda(NLFV lambda_) {lambda = lambda_;}

 private:
  arma::mat dnl;
  arma::mat dnli;
  NLFV lambda;
};


