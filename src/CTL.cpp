#include "CTL.h"

/*
 * Module for control options of optimization routines
*/
RCPP_MODULE(CTL){
  Rcpp::class_<CTRL>( "CTRL" )
    .constructor<int, double, double, double, bool>("sets the CTRL-options")

    .property("maxiters", &CTRL::get_maxiters, &CTRL::set_maxiters, "Getter and setter for maxiters")
    .property("abstol", &CTRL::get_abstol, &CTRL::set_abstol, "Getter and setter for abstol")
    .property("reltol", &CTRL::get_reltol, &CTRL::set_reltol, "Getter and setter for reltol")
    .property("feastol", &CTRL::get_feastol, &CTRL::set_feastol, "Getter and setter for feastol")
    .property("trace", &CTRL::get_trace, &CTRL::set_trace, "Getter and setter for trace")
    ;
}
