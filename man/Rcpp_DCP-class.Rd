\name{Rcpp_DCP-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_DCP-class}

\title{Class \code{"Rcpp_DCP"}}

\description{
  Class for definition of convex programs with non-linear constraints.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{x0}:}{Object of class \code{activeBindingFunction}:
      Initial values.}
    \item{\code{cList}:}{Object of class \code{activeBindingFunction}:
      Inequality constraints, class \code{CONEC}.}
    \item{\code{nList}:}{Object of class \code{activeBindingFunction}:
      List with elements of functions for evaluating non-linear
      constraints, their associated gradients and their associated
      Hessians.}
    \item{\code{A}:}{Object of class \code{activeBindingFunction}:
      Left-hand side of equality cosntraints.}
    \item{\code{b}:}{Object of class \code{activeBindingFunction}:
      Right-hand side of equality cosntraints.}
  }
}

\examples{
showClass("Rcpp_DCP")
}
\keyword{classes}
