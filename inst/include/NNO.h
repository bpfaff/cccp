/*
 * Class definition and methods for nonlinear constraints
*/
// Class definition for vectors/variables
class NNOV {
public:

  // constructors
  NNOV() : u(arma::mat()), dims(0L) {}
  NNOV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
  NNOV(arma::mat u_, int dims_): u( u_), dims( dims_) {}

  arma::mat u;
  int dims;
};
// Class definition for NT-scaling and Lagrange-Multipliers
class NNOS {
public:

  // constructors
  NNOS(arma::mat d_, arma::mat di_, NNOV lambda_): d(d_), di(di_), lambda(lambda_) {}

  arma::mat d;
  arma::mat di;
  NNOV lambda;
};
