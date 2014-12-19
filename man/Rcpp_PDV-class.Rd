\name{Rcpp_PDV-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_PDV-class}

\title{Class \code{"Rcpp_PDV"}}

\description{
  Class for primal/dual variables in convex programs.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{x}:}{Object of class \code{activeBindingFunction}:
      Primal variables.} 
    \item{\code{y}:}{Object of class \code{activeBindingFunction}: Dual
      variables.}
    \item{\code{s}:}{Object of class \code{activeBindingFunction}:
      Primal slack variables.}
    \item{\code{z}:}{Object of class \code{activeBindingFunction}: Dual
      slack variables.}
    \item{\code{kappa}:}{Object of class \code{activeBindingFunction}:
      Self-dual embedding variable; used in LPs, only.}
    \item{\code{tau}:}{Object of class \code{activeBindingFunction}:
      Self-dual embedding variable; used in LPs, only.}
  }
}

\examples{
showClass("Rcpp_PDV")
}
\keyword{classes}

