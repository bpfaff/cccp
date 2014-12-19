\name{Rcpp_CTRL-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_CTRL-class}

\title{Class \code{"Rcpp_CTRL"}}

\description{
  Class for control options used in optimization routines.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{maxiters}:}{Object of class
      \code{activeBindingFunction}: Integer, maximum count of iterations.}
    \item{\code{abstol}:}{Object of class \code{activeBindingFunction}:
      Numeric, absolute tolerance.}
    \item{\code{reltol}:}{Object of class \code{activeBindingFunction}:
      Numeric, relative tolerance.}
    \item{\code{feastol}:}{Object of class \code{activeBindingFunction}:
      Numeric, feasible tolerance.}
    \item{\code{trace}:}{Object of class \code{activeBindingFunction}:
      Logical, tracing progress during iterations.}
  }
}

\examples{
showClass("Rcpp_CTRL")
}
\keyword{classes}
