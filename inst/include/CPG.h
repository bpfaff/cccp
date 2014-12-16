#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(CTRL)
RCPP_EXPOSED_CLASS(CONEC)
RCPP_EXPOSED_CLASS(PDV)
RCPP_EXPOSED_CLASS(DQP)
RCPP_EXPOSED_CLASS(CPS)

#ifndef ARMA_H
#define ARMA_H
#include <RcppArmadillo.h>
#endif

/*
 * Class definition and methods for controlling optimization routines
*/
class CTRL {
 public:

  // constructors
 CTRL(int maxiters_, double abstol_, double reltol_, double feastol_, bool trace_): \
maxiters(maxiters_), abstol(abstol_), reltol(reltol_), feastol(feastol_), trace(trace_) {}
  // members
  int get_maxiters() {return maxiters;}
  void set_maxiters(int maxiters_) {maxiters = maxiters_;}
  double get_abstol() {return abstol;}
  void set_abstol(double abstol_) {abstol = abstol_;}
  double get_reltol() {return reltol;}
  void set_reltol(double reltol_) {reltol = reltol_;}
  double get_feastol() {return feastol;}
  void set_feastol(double feastol_) {feastol = feastol_;}
  bool get_trace() {return trace;}
  void set_trace(bool trace_) {trace = trace_;}

 private:
  int maxiters;
  double abstol;
  double reltol;
  double feastol;
  bool trace;
};

/*
 * Class definition for inequality (cone) constraints
*/
class CONEC {
 public:
 CONEC() : conTypes(std::vector<std::string>()), Gmats(std::vector<arma::mat>()), \
    hvecs(std::vector<arma::vec>()), dims(arma::uvec()), K(0) {}
 CONEC(Rcpp::CharacterVector conTypes_, Rcpp::List Gmats_, Rcpp::List hvecs_, Rcpp::IntegerVector dims_, int K_): \
  conTypes(), Gmats(), hvecs(), K(K_){
    conTypes = Rcpp::as< std::vector<std::string> >(conTypes_);
    Gmats = Rcpp::as< std::vector<arma::mat> >(Gmats_);
    hvecs = Rcpp::as< std::vector<arma::vec> >(hvecs_);
    dims = Rcpp::as<arma::uvec>(dims_);
  }
  // members
  std::vector<std::string> get_conTypes() {return conTypes;}
  void set_conTypes(std::vector<std::string> conTypes_) {conTypes = conTypes_;}
  std::vector<arma::mat> get_Gmats() {return Gmats;}
  void set_Gmats(std::vector<arma::mat> Gmats_) {Gmats = Gmats_;}
  std::vector<arma::vec> get_hvecs() {return hvecs;}
  void set_hvecs(std::vector<arma::vec> hvecs_) {hvecs = hvecs_;}
  arma::uvec get_dims() {return dims;}
  void set_dims(arma::uvec dims_) {dims = dims_;}
  int get_K() {return K;}
  void set_K(int K_) {K = K_;}

  friend class DQP;

 private:
  std::vector<std::string> conTypes;
  std::vector<arma::mat> Gmats;
  std::vector<arma::vec> hvecs;
  arma::uvec dims;
  int K;
};

/*
 * Class for definition of Quadratic program
*/
class DQP {
 public:

  // constructors
 DQP() : P(arma::mat()), q(arma::vec()), A(arma::mat()), b(arma::vec()), cList(CONEC()) {}
  DQP(arma::mat P_, arma::vec q_, arma::mat A_, arma::vec b_, CONEC cList_): \
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
  CONEC get_cList() {return cList;}
  void set_cList(CONEC cList_) {cList = cList_;}

  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  arma::mat rprim(PDV& pdv);
  std::vector<arma::mat> rcent(PDV& pdv);
  arma::mat rdual(PDV& pdv);
  CPS* cps(CTRL& ctrl);

 private:
  arma::mat P;
  arma::vec q;
  arma::mat A;
  arma::vec b;
  CONEC cList;
};

/*
 * Class definition for primal/dual variables
*/
class PDV {
 public:

  // constructors
 PDV() : x(arma::mat()), y(arma::mat()), s(std::vector<arma::mat>()), z(std::vector<arma::mat>()), \
    kappa(1.0), tau(1.0) {}
 PDV(arma::mat x_, arma::mat y_, std::vector<arma::mat> s_, std::vector<arma::mat> z_, \
     double kappa_, double tau_):  \
  x(x_), y(y_), s(s_), z(z_), kappa(kappa_), tau(tau_) {}
  // members
  arma::mat get_x() {return x;}
  void set_x(arma::mat x_) {x = x_;}
  arma::mat get_y() {return y;}
  void set_y(arma::mat y_) {y = y_;}
  std::vector<arma::mat> get_s() {return s;}
  void set_s(std::vector<arma::mat> s_) {s = s_;}
  std::vector<arma::mat> get_z() {return z;}
  void set_z(std::vector<arma::mat> z_) {z = z_;}
  double get_kappa() {return kappa;}
  void set_kappa(double kappa_) {kappa = kappa_;}
  double get_tau() {return tau;}
  void set_tau(double tau_) {tau = tau_;}

  friend class DQP;

 private:
  arma::mat x;
  arma::mat y;
  std::vector<arma::mat> s;
  std::vector<arma::mat> z;
  double kappa;
  double tau;
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

