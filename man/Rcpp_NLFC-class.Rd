\name{Rcpp_NLFC-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_NLFC-class}

\title{Class \code{"Rcpp_NLFC"}}

\description{
  Class for nonlinear constraints. 
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{G}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{h}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{dims}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
  }
}

\examples{
showClass("Rcpp_NLFC")
}
\keyword{classes}
