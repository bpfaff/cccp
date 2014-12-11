\name{Rcpp_NLFV-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_NLFV-class}

\title{Class \code{"Rcpp_NLFV"}}

\description{
  Class for variables that are related to nonlinear constraints.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{dims}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{u}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
    \item{\code{ntsc(...)}:}{ ~~ }
    \item{\code{uone(...)}:}{ ~~ }
    \item{\code{udot(...)}:}{ ~~ }
    \item{\code{uprd(...)}:}{ ~~ }
    \item{\code{uinv(...)}:}{ ~~ }
    \item{\code{umss(...)}:}{ ~~ }
    \item{\code{umsa(...)}:}{ ~~ }
  }
}

\examples{
showClass("Rcpp_NLFV")
}
\keyword{classes}
