#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(CTRL)
RCPP_EXPOSED_CLASS(CONEC)
RCPP_EXPOSED_CLASS(PDV)
RCPP_EXPOSED_CLASS(DQP)
RCPP_EXPOSED_CLASS(DLP)
RCPP_EXPOSED_CLASS(DNL)
RCPP_EXPOSED_CLASS(DCP)
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
 CTRL(): params(Rcpp::List::create()) {}
 CTRL(Rcpp::List params_): params(params_){}
  // members
  Rcpp::List get_params() {return params;}
  void set_params(Rcpp::List params_) {params = params_;}
  Rcpp::List params;
};

/*
 * Class definition for inequality (cone) constraints
*/
class CONEC {
 public:
 CONEC() : cone(std::vector<std::string>()), G(mat()), 
    h(mat()), sidx(umat()), dims(uvec()), K(0), n(0) {}
 CONEC(std::vector<std::string> cone_, mat G_, mat h_, umat sidx_, uvec dims_, int K_, int n_): 
  cone(cone_), G(G_), h(h_), sidx(sidx_), dims(dims_), K(K_), n(n_){}
 CONEC(int n_): cone(std::vector<std::string>()), G(mat()), 
    h(mat()), sidx(umat()), dims(uvec()), K(0), n(n_){}
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
  int get_n() {return n;}
  void set_n(int n_) {n = n_;}

  friend class DCP;
  friend class DLP;
  friend class DNL;
  friend class DQP;
  friend CPS* rpp(mat x0, mat P, mat mrc, CTRL& ctrl);
  friend CPS* gpp(std::vector<mat> FList, std::vector<mat> gList, CONEC& cList, mat A, mat b, CTRL& ctrl);

  double snrm2(mat s);
  vec sdot(mat s, mat z);
  vec smss(mat u);
  mat sone();
  mat sprd(mat s, mat z);
  mat sinv(mat s, mat z);
  mat sams1(mat u, double alpha);
  mat sslb(mat s, mat lambda, bool invers);
  mat ssnt(mat s, std::vector<std::map<std::string,mat> > WList, 
	   bool invers, bool transp);
  mat getLambda(std::vector<std::map<std::string,mat> > WList);
  mat gwwg(std::vector<std::map<std::string,mat> > WList);
  mat gwwz(std::vector<std::map<std::string,mat> > WList, mat z);
  mat SorZupdate(mat SorZ, mat Lambda, double step);
  PDV* initpdv(int p);
  std::vector<std::map<std::string,mat> > initnts();
  std::vector<std::map<std::string,mat> > ntsc(mat s, mat z);
  std::vector<std::map<std::string,mat> > 
    ntsu(mat s, mat z, 
	 std::vector<std::map<std::string,mat> > WList);

 private:
  std::vector<std::string> cone;
  mat G;
  mat h;
  umat sidx;
  uvec dims;
  int K;
  int n;
};
/*
 * Class for definition of Quadratic programs
*/
class DQP {
 public:

  // constructors
 DQP() : P(mat()), q(vec()), A(mat()), b(vec()), cList(CONEC()) {}
  DQP(mat P_, vec q_, mat A_, vec b_, CONEC cList_): 
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

  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  mat rprim(PDV& pdv);
  mat rcent(PDV& pdv);
  mat rdual(PDV& pdv);
  PDV* sxyz(PDV* pdv, mat LHS, mat RHS, 
	    std::vector<std::map<std::string,mat> > WList);
  CPS* cps(CTRL& ctrl);

 private:
  mat P;
  vec q;
  mat A;
  vec b;
  CONEC cList;
};

/*
 * Class for definition of Linear programs
*/
class DLP {
 public:

  // constructors
 DLP() : q(vec()), A(mat()), b(vec()), cList(CONEC()) {}
  DLP(vec q_, mat A_, vec b_, CONEC cList_): 
  q(q_), A(A_), b(b_), cList(cList_) {}
  // members
  vec get_q() {return q;}
  void set_q(vec q_) {q = q_;}
  mat get_A() {return A;}
  void set_A(mat A_) {A = A_;}
  vec get_b() {return b;}
  void set_b(vec b_) {b = b_;}
  CONEC get_cList() {return cList;}
  void set_cList(CONEC cList_) {cList = cList_;}

  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  mat rprim(PDV& pdv);
  mat rcent(PDV& pdv);
  mat rdual(PDV& pdv);
  PDV* sxyz(PDV* pdv, mat LHS, mat RHS, 
	    std::vector<std::map<std::string,mat> > WList);
  CPS* cps(CTRL& ctrl);

 private:
  vec q;
  mat A;
  vec b;
  CONEC cList;
};

/*
 * Class for definition of Linear programs with non-linear constraints
*/
class DNL {
 public:

  // constructors
 DNL() : q(vec()), A(mat()), b(vec()), cList(CONEC()), x0(mat()), 
    nList(Rcpp::List::create()) {}
 DNL(vec q_, mat A_, vec b_, CONEC cList_, mat x0_, Rcpp::List nList_): 
  q(q_), A(A_), b(b_), cList(cList_), x0(x0_), nList(nList_) {}
  // members
  vec get_q() {return q;}
  void set_q(vec q_) {q = q_;}
  mat get_A() {return A;}
  void set_A(mat A_) {A = A_;}
  vec get_b() {return b;}
  void set_b(vec b_) {b = b_;}
  CONEC get_cList() {return cList;}
  void set_cList(CONEC cList_) {cList = cList_;}
  mat get_x0() {return x0;}
  void set_x0(mat x0_) {x0 = x0_;}
  Rcpp::List get_nList() {return nList;}
  void set_nList(Rcpp::List nList_) {nList = nList_;}

  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  mat rprim(PDV& pdv);
  mat rcent(PDV& pdv);
  mat rdual(PDV& pdv);
  PDV* sxyz(PDV* pdv, mat LHS, mat RHS, 
	    std::vector<std::map<std::string,mat> > WList);
  CPS* cps(CTRL& ctrl);

 private:
  vec q;
  mat A;
  vec b;
  CONEC cList;
  mat x0;
  Rcpp::List nList;
};


/*
 * Class for definition of convex programs with non-linear constraints
*/
class DCP {
 public:

  // constructors
 DCP() : x0(mat()), cList(CONEC()), nList(Rcpp::List::create()), 
    A(mat()), b(vec()) {}
 DCP(mat x0_, CONEC cList_, Rcpp::List nList_, mat A_, vec b_): 
  x0(x0_), cList(cList_), nList(nList_), A(A_), b(b_) {}
  // members
  mat get_x0() {return x0;}
  void set_x0(mat x0_) {x0 = x0_;}
  CONEC get_cList() {return cList;}
  void set_cList(CONEC cList_) {cList = cList_;}
  Rcpp::List get_nList() {return nList;}
  void set_nList(Rcpp::List nList_) {nList = nList_;}
  mat get_A() {return A;}
  void set_A(mat A_) {A = A_;}
  vec get_b() {return b;}
  void set_b(vec b_) {b = b_;}

  double pobj(PDV& pdv);
  double dobj(PDV& pdv);
  double certp(PDV& pdv);
  double certd(PDV& pdv);
  mat rprim(PDV& pdv);
  mat rcent(PDV& pdv);
  mat rdual(PDV& pdv);
  PDV* sxyz(PDV* pdv, mat LHS, 
	    std::vector<std::map<std::string,mat> > WList);
  CPS* cps(CTRL& ctrl);

 private:
  mat x0;
  CONEC cList;
  Rcpp::List nList;
  mat A;
  vec b;
};


/*
 * Class definition for primal/dual variables
*/
class PDV {
 public:

  // constructors
 PDV() : x(mat()), y(mat()), s(mat()), z(mat()), kappa(1.0), tau(1.0) {}
 PDV(mat x_, mat y_, mat s_, mat z_, double kappa_, double tau_):  
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

  friend class DCP;
  friend class DLP;
  friend class DNL;
  friend class DQP;
  friend class CONEC;
  friend CPS* rpp(mat x0, mat P, mat mrc, CTRL& ctrl);
  friend CPS* gpp(std::vector<mat> FList, std::vector<mat> gList, CONEC& cList, mat A, mat b, CTRL& ctrl);

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
 CPS() : pdv(PDV()), state(Rcpp::NumericVector::create()), 
    status("unknown"), niter(0), sidx(umat()) 
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
 CPS(PDV pdv_, Rcpp::NumericVector state_, Rcpp::String status_, 
     int niter_, umat sidx_): 
  pdv(pdv_), state(state_), status(status_), niter(niter_) , sidx(sidx_){}
  // members
  PDV get_pdv() {return pdv;}
  void set_pdv(PDV pdv_) {pdv = pdv_;}
  Rcpp::NumericVector get_state() {return state;}
  void set_state(Rcpp::NumericVector state_) {state = state_;}
  Rcpp::String get_status() {return status;}
  void set_status(Rcpp::String status_) {status = status_;}
  int get_niter() {return niter;}
  void set_niter(int niter_) {niter = niter_;}
  umat get_sidx() {return sidx;}
  void set_sidx(umat sidx_) {sidx = sidx_;}

  PDV pdv;

 private:
  Rcpp::NumericVector state;
  Rcpp::String status;
  int niter;
  umat sidx;
};

// Function for solving risk parity portfolios
CPS* rpp(mat x0, mat P, mat mrc, CTRL& ctrl);
// Function for solving geometric programs
CPS* gpp(std::vector<mat> FList, std::vector<mat> gList, CONEC& cList, mat A, mat b, CTRL& ctrl);
