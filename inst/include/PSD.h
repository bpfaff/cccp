/*
 * Class definition and methods for positive semi-definite cone constraints
*/
class PSDV {
public:

  // constructors
  PSDV() : u(arma::mat()), dims(0L) {}
  PSDV(int dims_) : u(arma::mat().ones(dims_, 1)), dims(dims_) {}
  PSDV(arma::mat u_, int dims_): u(u_), dims(dims_) {}

  arma::mat u;	
  int dims;
};
