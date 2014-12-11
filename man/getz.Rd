\name{getz}
\alias{getz}
\alias{getz,Rcpp_PDV-method}
\alias{getz,Rcpp_CPS-method}

\title{
  Extractor for slot \code{s}
}

\description{
  Returns the value of the list-variable(s) \code{z}.
}

\usage{
\S4method{getz}{Rcpp_PDV}(object)  
\S4method{getz}{Rcpp_CPS}(object)  
}

\arguments{
  \item{object}{An object of either reference-class \code{Rcpp_PDV} or
  \code{Rcpp_CPS}.} 
}

\value{
  \code{list}, the value of the list-variable(s) \code{z} (dual-slack).
}

\keyword{optimize}
