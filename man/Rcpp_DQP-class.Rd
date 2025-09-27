\name{Rcpp_DQP-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_DQP-class}

\title{Class \code{"Rcpp_DQP"}}

\description{
  Class for definition of quadratic programs.
}

\section{Extends}{
Class \link[Rcpp]{C++Object-class}, directly.
All reference classes extend and inherit methods from \link{envRefClass-class}.
}
       
\section{Fields}{
  \describe{
    \item{\code{P}:}{Object of class \code{activeBindingFunction}:
      Matrix related to quadratic term in objective.}
    \item{\code{q}:}{Object of class \code{activeBindingFunction}:
      Matrix related to linear term in objective.}
    \item{\code{A}:}{Object of class \code{activeBindingFunction}:
      Left-hand side of equality cosntraints.}
    \item{\code{b}:}{Object of class \code{activeBindingFunction}:
      Right-hand side of equality cosntraints.}
    \item{\code{cList}:}{Object of class \code{activeBindingFunction}:
      Inequality constraints, class \code{CONEC}.}
  }
}

\examples{
showClass("Rcpp_DQP")
}
\keyword{classes}
