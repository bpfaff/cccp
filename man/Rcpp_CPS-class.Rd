\name{Rcpp_CPS-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_CPS-class}

\title{Class \code{"Rcpp_CPS"}}

\description{
  Class for solution of convex programs.
}

\section{Extends}{
Class \link[Rcpp]{C++Object-class}, directly.
All reference classes extend and inherit methods from \link{envRefClass-class}.
}
       
\section{Fields}{
  \describe{
    \item{\code{pdv}:}{Object of class \code{activeBindingFunction}:
      Primal-dual variables.}
    \item{\code{state}:}{Object of class \code{activeBindingFunction}:
      Vector of state variables in convex programs.}
    \item{\code{status}:}{Object of class \code{activeBindingFunction}:
      Character indicating the status of the returned solution.}
    \item{\code{niter}:}{Object of class \code{activeBindingFunction}:
      Integer, count of iterations.}
    \item{\code{sidx}:}{Object of class \code{activeBindingFunction}:
      Integer matrix, start and end indices of slack variables.}
  }
}

\examples{
showClass("Rcpp_CPS")
}
\keyword{classes}
