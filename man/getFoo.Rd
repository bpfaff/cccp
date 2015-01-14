\name{getFoo}
\alias{gets}
\alias{gets,Rcpp_PDV-method}
\alias{gets,Rcpp_CPS-method}
\alias{getx}
\alias{getx,Rcpp_PDV-method}
\alias{getx,Rcpp_CPS-method}
\alias{gety}
\alias{gety,Rcpp_PDV-method}
\alias{gety,Rcpp_CPS-method}
\alias{getz}
\alias{getz,Rcpp_PDV-method}
\alias{getz,Rcpp_CPS-method}
\alias{getstate}
\alias{getstate,Rcpp_CPS-method}
\alias{getstatus}
\alias{getstatus,Rcpp_CPS-method}
\alias{getniter}
\alias{getniter,Rcpp_CPS-method}
\alias{getparams}
\alias{getparams,Rcpp_CTRL-method}

\title{
  Extractor methods for reference class objects
}

\description{
  Returns a member of reference class objects.
}

\usage{
\S4method{getx}{Rcpp_PDV}(object)  
\S4method{getx}{Rcpp_CPS}(object)  
\S4method{gety}{Rcpp_PDV}(object)  
\S4method{gety}{Rcpp_CPS}(object)  
\S4method{gets}{Rcpp_PDV}(object)  
\S4method{gets}{Rcpp_CPS}(object)  
\S4method{getz}{Rcpp_PDV}(object)  
\S4method{getz}{Rcpp_CPS}(object)  
\S4method{getstate}{Rcpp_CPS}(object)  
\S4method{getstatus}{Rcpp_CPS}(object)  
\S4method{getniter}{Rcpp_CPS}(object)  
\S4method{getparams}{Rcpp_CTRL}(object)  
}

\arguments{
  \item{object}{An object of either reference-class \code{Rcpp_PDV} or
  \code{Rcpp_CPS}, or \code{Rcpp_CTRL}.} 
}

\value{
  The relevant member object of the class.
}

\keyword{optimize}
