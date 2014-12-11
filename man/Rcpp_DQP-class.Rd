\name{Rcpp_DQP-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_DQP-class}
\alias{show}

\title{Class \code{"Rcpp_DQP"}}

\description{
  Class for definition of quadratic programs.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
       
\section{Fields}{
  \describe{
    \item{\code{P}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{q}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{A}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{b}:}{Object of class \code{activeBindingFunction} ~~ }
    \item{\code{cList}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
    \item{\code{pobj()}:}{ ~~ }
    \item{\code{cps()}:}{ ~~ }
  }
}


\examples{
showClass("Rcpp_DQP")
}
\keyword{classes}
