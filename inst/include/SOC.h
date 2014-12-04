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

  arma::mat u;	
  int dims;
};
// Class definition for NT-scaling and Lagrange-Multipliers
class SOCS {
public:

  // constructors
  SOCS(arma::mat v_, double beta_, SOCV lambda_): v(v_), beta(beta_), lambda(lambda_) {}

  arma::mat v;
  double beta;
  SOCV lambda;
};
