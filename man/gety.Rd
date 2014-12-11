\name{gety}
\alias{gety}
\alias{gety,Rcpp_PDV-method}
\alias{gety,Rcpp_CPS-method}

\title{
  Extractor for slot \code{y}
}

\description{
  Returns the value of the dual variable \code{y}.
}

\usage{
\S4method{gety}{Rcpp_PDV}(object)  
\S4method{gety}{Rcpp_CPS}(object)  
}

\arguments{
  \item{object}{An object of either reference-class \code{Rcpp_PDV} or
  \code{Rcpp_CPS}.} 
}

\value{
  \code{numeric}, the value of the dual-variable(s) \code{y}.
}

\keyword{optimize}
