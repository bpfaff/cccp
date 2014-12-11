#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(PDV)
RCPP_EXPOSED_CLASS(DQP)

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
 PDV(arma::vec x_, arma::vec y_, Rcpp::List s_, Rcpp::List z_, double kappa_, double tau_):  \
   x(x_), y(y_), s(s_), z(z_), kappa(kappa_), tau(tau_) {}
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


/*
 * Class for definition of Quadratic program
*/
class DQP {
 public:

  // constructors
 DQP() : P(arma::mat()), q(arma::vec()), A(arma::mat()), b(arma::vec()), cList(Rcpp::List::create()) {}
 DQP(arma::mat P_, arma::vec q_, arma::mat A_, arma::vec b_, Rcpp::List cList_):  \
   P(P_), q(q_), A(A_), b(b_), cList(cList_) {}
   // members
   arma::mat get_P() {return P;}
   void set_P(arma::mat P_) {P = P_;}
   arma::vec get_q() {return q;}
   void set_q(arma::vec q_) {q = q_;}
   arma::mat get_A() {return A;}
   void set_A(arma::mat A_) {A = A_;}
   arma::vec get_b() {return b;}
   void set_b(arma::vec b_) {b = b_;}
   Rcpp::List get_cList() {return cList;}
   void set_cList(Rcpp::List cList_) {cList = cList_;}

 private:
   arma::mat P;
   arma::vec q;
   arma::mat A;
   arma::vec b;
   Rcpp::List cList;
};

/*
 * Class for solution of convex programs
*/
class CPS {
 public:

  // constructors
 CPS() : pdv(PDV()), state(Rcpp::NumericVector::create()), status("unknown"), niter(0) 
    {
      state["pobj"] = NA_REAL;
      state["dobj"] = NA_REAL;
      state["dgap"] = NA_REAL;
      state["rdgap"] = NA_REAL;
      state["certp"] = NA_REAL;
      state["certd"] = NA_REAL;
      state["pslack"] = NA_REAL;
      state["dslack"] = NA_REAL;
      status = "unknown";
    }
 CPS(PDV pdv_, Rcpp::NumericVector state_, Rcpp::String status_, int niter_):	\
  pdv(pdv_), state(state_), status(status_), niter(niter_) {}
  // members
  PDV get_pdv() {return pdv;}
  void set_pdv(PDV pdv_) {pdv = pdv_;}
  Rcpp::NumericVector get_state() {return state;}
  void set_state(Rcpp::NumericVector state_) {state = state_;}
  Rcpp::String get_status() {return status;}
  void set_status(Rcpp::String status_) {status = status_;}
  int get_niter() {return niter;}
  void set_niter(int niter_) {niter = niter_;}

 private:
  PDV pdv;
  Rcpp::NumericVector state;
  Rcpp::String status;
  int niter;
};

