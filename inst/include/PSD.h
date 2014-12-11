#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(PSDV)
RCPP_EXPOSED_CLASS(PSDS)
RCPP_EXPOSED_CLASS(PSDC)

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
 PSDV(int dims_) : u(arma::mat(dims_, dims_)), dims(dims_) 
  {
    u.eye();
    u.reshape(dims * dims, 1);
  }
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
  PSDS* ntsc(const PSDV& z) const;

 private:
  arma::mat u;	
  int dims;
};

// Class definition for NT-scaling and Lagrange-Multipliers
class PSDS {
 public:

  // constructors
 PSDS(arma::mat r_, arma::mat rti_, PSDV lambda_): r(r_), rti(rti_), lambda(lambda_) {}

  arma::mat get_r() {return r;}
  void set_r(arma::mat r_) {r = r_;}
  arma::mat get_rti() {return rti;}
  void set_rti(arma::mat rti_) {rti = rti_;}
  PSDV get_lambda() {return lambda;}
  void set_lambda(PSDV lambda_) {lambda = lambda_;}

 private:
  arma::mat r;
  arma::mat rti;
  PSDV lambda;
};


// Class definition for positive semidefinite constraints
class PSDC {
 public:

  // constructors
 PSDC(arma::mat G_, PSDV h_, int dims_): G(G_), h(h_), dims(dims_) {}
  // members
  arma::mat get_G() {return G;}
  void set_G(arma::mat G_) {G = G_;}
  PSDV get_h() {return h;}
  void set_h(PSDV h_) {h = h_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

 private:
  arma::mat G;
  PSDV h;
  int dims;
};
