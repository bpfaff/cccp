\name{gets}
\alias{gets}
\alias{gets,Rcpp_PDV-method}
\alias{gets,Rcpp_CPS-method}

\title{
  Extractor for slot \code{s}
}

\description{
  Returns the value of the list-variable(s) \code{s}.
}

\usage{
\S4method{gets}{Rcpp_PDV}(object)  
\S4method{gets}{Rcpp_CPS}(object)  
}

\arguments{
  \item{object}{An object of either reference-class \code{Rcpp_PDV} or
  \code{Rcpp_CPS}.} 
}

\value{
  \code{list}, the value of the list-variable(s) \code{s} (primal slack).
}

\keyword{optimize}
