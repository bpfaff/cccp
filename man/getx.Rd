\name{getx}
\alias{getx}
\alias{getx,Rcpp_PDV-method}
\alias{getx,Rcpp_CPS-method}

\title{
  Extractor for slot \code{x}
}

\description{
  Returns the value of the objective variable \code{x}.
}

\usage{
\S4method{getx}{Rcpp_PDV}(object)  
\S4method{getx}{Rcpp_CPS}(object)  
}

\arguments{
  \item{object}{An object of either reference-class \code{Rcpp_PDV} or
  \code{Rcpp_CPS}.} 
}

\value{
  \code{numeric}, the value of the objective-variable(s) \code{x}.
}

\keyword{optimize}
