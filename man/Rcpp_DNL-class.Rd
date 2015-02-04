\name{Rcpp_DNL-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_DNL-class}

\title{Class \code{"Rcpp_DNL"}}

\description{
  Class for definition of linear programs with non-linear constraints.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{q}:}{Object of class \code{activeBindingFunction}:
      Matrix related to linear term in objective.}
    \item{\code{A}:}{Object of class \code{activeBindingFunction}:
      Left-hand side of equality cosntraints.}
    \item{\code{b}:}{Object of class \code{activeBindingFunction}:
      Right-hand side of equality cosntraints.}
    \item{\code{cList}:}{Object of class \code{activeBindingFunction}:
      Inequality constraints, class \code{CONEC}.}
    \item{\code{x0}:}{Object of class \code{activeBindingFunction}:
      Initial values.}
    \item{\code{nList}:}{Object of class \code{activeBindingFunction}:
      List with elements of functions for evaluating non-linear
      constraints, their associated gradients and their associated
      Hessians.}
  }
}

\examples{
showClass("Rcpp_DNL")
}
\keyword{classes}
