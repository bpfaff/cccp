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
    \item{\code{cone}:}{Object of class \code{activeBindingFunction}:
      Type of cone constraints.}
    \item{\code{G}:}{Object of class \code{activeBindingFunction}:
      Left-hand side of inequality constraints.}
    \item{\code{h}:}{Object of class \code{activeBindingFunction}:
      Right-hand side of inequality constraints.}
    \item{\code{sidx}:}{Object of class \code{activeBindingFunction}:
      Row index for subsets of cone constraints.}
    \item{\code{dims}:}{Object of class \code{activeBindingFunction}:
      Dimension of cone constraints.}  
    \item{\code{K}:}{Object of class \code{activeBindingFunction}: Count
      of inequality constraints.}  
    \item{\code{n}:}{Object of class \code{activeBindingFunction}: Count
      of variables in objective.}  
  }
}

\examples{
showClass("Rcpp_CONEC")
}
\keyword{classes}
