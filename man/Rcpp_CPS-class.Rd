\name{Rcpp_CPS-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_CPS-class}

\title{Class \code{"Rcpp_CPS"}}

\description{
  Class for solution of convex programs.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{pdv}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{state}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{status}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{niter}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}

\examples{
showClass("Rcpp_CPS")
}
\keyword{classes}
