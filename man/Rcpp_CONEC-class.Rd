\name{Rcpp_CONEC-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_CONEC-class}

\title{Class \code{"Rcpp_CONEC"}}

\description{
  Class for inequality (cone) constraints.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{coneTypes}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{Gmats}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{hvecs}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{K}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}

\examples{
showClass("Rcpp_CONEC")
}
\keyword{classes}
