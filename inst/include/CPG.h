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
using namespace arma;

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
 CONEC() : cone(std::vector<std::string>()), G(mat()), \
    h(mat()), sidx(umat()), dims(uvec()), K(0) {}
 CONEC(std::vector<std::string> cone_, mat G_, mat h_, umat sidx_, uvec dims_, int K_): \
  cone(cone_), G(G_), h(h_), sidx(sidx_), dims(dims_), K(K_){}
  // members
  std::vector<std::string> get_cone() {return cone;}
  void set_cone(std::vector<std::string> cone_) {cone = cone_;}
  mat get_G() {return G;}
  void set_G(mat G_) {G = G_;}
  mat get_h() {return h;}
  void set_h(mat h_) {h = h_;}
  umat get_sidx() {return sidx;}
  void set_sidx(umat sidx_) {sidx = sidx_;}
  uvec get_dims() {return dims;}
  void set_dims(uvec dims_) {dims = dims_;}
  int get_K() {return K;}
  void set_K(int K_) {K = K_;}

  friend class DQP;

 private:
  std::vector<std::string> cone;
  mat G;
  mat h;
  umat sidx;
  uvec dims;
  int K;
};
/*
 * Class for definition of Quadratic program
*/
class DQP {
 public:

  // constructors
 DQP() : P(mat()), q(vec()), A(mat()), b(vec()), cList(CONEC()) {}
  DQP(mat P_, vec q_, mat A_, vec b_, CONEC cList_): \
  P(P_), q(q_), A(A_), b(b_), cList(cList_) {}
  // members
  mat get_P() {return P;}
  void set_P(mat P_) {P = P_;}
  vec get_q() {return q;}
  void set_q(vec q_) {q = q_;}
  mat get_A() {return A;}
  void set_A(mat A_) {A = A_;}
  vec get_b() {return b;}
  void set_b(vec b_) {b = b_;}
  CONEC get_cList() {return cList;}
  void set_cList(CONEC cList_) {cList = cList_;}

  double snrm2(mat s);
  vec sdot(mat s, mat z);
  vec smss(mat s);
  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  mat rprim(PDV& pdv);
  mat rcent(PDV& pdv);
  mat rdual(PDV& pdv);
  mat sams1(mat s, double alpha);
  PDV* initpdv();
  std::vector<std::map<std::string,mat> > initnts();
  mat gwwg(std::vector<std::map<std::string,mat> > WList);
  mat gwwz(std::vector<std::map<std::string,mat> > WList, mat z);
  PDV* sxyz(PDV* pdv, mat LHS, mat RHS, std::vector<std::map<std::string,mat> > WList);
  mat ssnt(mat s, std::vector<std::map<std::string,mat> > WList, bool invers, bool transp);
  CPS* cps(CTRL& ctrl);

 private:
  mat P;
  vec q;
  mat A;
  vec b;
  CONEC cList;
};

/*
 * Class definition for primal/dual variables
*/
class PDV {
 public:

  // constructors
 PDV() : x(mat()), y(mat()), s(mat()), z(mat()), kappa(1.0), tau(1.0) {}
 PDV(mat x_, mat y_, mat s_, mat z_, double kappa_, double tau_):  \
  x(x_), y(y_), s(s_), z(z_), kappa(kappa_), tau(tau_) {}

  // members
  mat get_x() {return x;}
  void set_x(mat x_) {x = x_;}
  mat get_y() {return y;}
  void set_y(mat y_) {y = y_;}
  mat get_s() {return s;}
  void set_s(mat s_) {s = s_;}
  mat get_z() {return z;}
  void set_z(mat z_) {z = z_;}
  double get_kappa() {return kappa;}
  void set_kappa(double kappa_) {kappa = kappa_;}
  double get_tau() {return tau;}
  void set_tau(double tau_) {tau = tau_;}

  friend class DQP;
  friend class CONEC;

 private:
  mat x;
  mat y;
  mat s;
  mat z;
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

