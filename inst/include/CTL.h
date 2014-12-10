#ifndef REFC_H
#define REFC_H
#include <RcppCommon.h>
#endif
// forward declarations and helping module classes 
RCPP_EXPOSED_CLASS(CTRL)

#ifndef RCPP_H
#define RCPP_H
#include <Rcpp.h>
#endif

/*
 * Class definition and methods for controlling optimization routines
*/
// Class definition for vectors/variables 
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


