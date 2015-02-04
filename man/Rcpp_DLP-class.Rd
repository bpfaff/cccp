\name{Rcpp_DLP-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_DLP-class}

\title{Class \code{"Rcpp_DLP"}}

\description{
  Class for definition of linear programs.
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
  }
}

\examples{
showClass("Rcpp_DLP")
}
\keyword{classes}
