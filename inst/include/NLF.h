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

  arma::mat u;
  int dims;
};
// Class definition for NT-scaling and Lagrange-Multipliers
class NLFS {
public:

  // constructors
  NLFS(arma::mat dnl_, arma::mat dnli_, NLFV lambda_): dnl(dnl_), dnli(dnli_), lambda(lambda_) {}

  arma::mat dnl;
  arma::mat dnli;
  NLFV lambda;
};
