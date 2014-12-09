#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(PSDV)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition and methods for positive semi-definite cone constraints
*/
class PSDV {
 public:

  // constructors
 PSDV() : u(arma::mat()), dims(0L) {}
 PSDV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
 PSDV(arma::mat u_, int dims_): u(u_), dims(dims_) {}
  // members
  arma::mat get_u() {return u;}
  void set_u(arma::mat u_) {u = u_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

  double umss() const;
  double udot(const PSDV& z) const;
  PSDV* uone() const;
  PSDV* uprd(const PSDV& z) const;
  PSDV* uinv(const PSDV& z) const;
  PSDV* umsa1(double alpha) const;
  PSDV* umsa2(double alpha, const PSDV& lambda, arma::vec sigma) const;

 private:
  arma::mat u;	
  int dims;
};
