#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(NNOV)
RCPP_EXPOSED_CLASS(NNOS)
RCPP_EXPOSED_CLASS(NNOC)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition and methods for nonlinear constraints
*/
// Class definition for vectors/variables 
class NNOV {
 public:

  // constructors
 NNOV() : u(arma::mat()), dims(0L) {}
 NNOV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
 NNOV(arma::mat u_, int dims_): u(u_), dims(dims_) {}
  // members
  arma::mat get_u() {return u;}
  void set_u(arma::mat u_) {u = u_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

  NNOV* uone() const;
  double umss() const;
  double udot(const NNOV& z) const;
  NNOV* uprd(const NNOV& z) const;
  NNOV* uinv(const NNOV& z) const;
  NNOV* umsa(double alpha, bool init) const;
  NNOS* ntsc(const NNOV& z) const;

 private:
  arma::mat u;
  int dims;
};

// Class definition for NT-scaling and Lagrange-Multipliers
class NNOS {
 public:

  // constructors
 NNOS(arma::mat d_, arma::mat di_, NNOV lambda_): d(d_), di(di_), lambda(lambda_) {}

  arma::mat get_d() {return d;}
  void set_d(arma::mat d_) {d = d_;}
  arma::mat get_di() {return di;}
  void set_di(arma::mat di_) {di = di_;}
  NNOV get_lambda() {return lambda;}
  void set_lambda(NNOV lambda_) {lambda = lambda_;}

 private:
  arma::mat d;
  arma::mat di;
  NNOV lambda;
};


// Class definition for linear constraints
class NNOC {
 public:

  // constructors
 NNOC(arma::mat G_, NNOV h_, int dims_): G(G_), h(h_), dims(dims_) {}
  // members
  arma::mat get_G() {return G;}
  void set_G(arma::mat G_) {G = G_;}
  NNOV get_h() {return h;}
  void set_h(NNOV h_) {h = h_;}
  int get_dims() {return dims;}
  void set_dims(int dims_) {dims = dims_;}

 private:
  arma::mat G;
  NNOV h;
  int dims;
};
