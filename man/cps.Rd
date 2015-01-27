\name{cps}
\alias{cps}
\alias{cps,Rcpp_DLP,Rcpp_CTRL-method}
\alias{cps,Rcpp_DNL,Rcpp_CTRL-method}
\alias{cps,Rcpp_DQP,Rcpp_CTRL-method}
\alias{cps,Rcpp_DCP,Rcpp_CTRL-method}

\title{
  Solving a convex program
}

\description{
  This function returns an optimal point for a cone constraint convex
  program. 
}

\usage{
\S4method{cps}{Rcpp_DCP,Rcpp_CTRL}(cpd, ctrl)  
\S4method{cps}{Rcpp_DLP,Rcpp_CTRL}(cpd, ctrl)  
\S4method{cps}{Rcpp_DNL,Rcpp_CTRL}(cpd, ctrl)  
\S4method{cps}{Rcpp_DQP,Rcpp_CTRL}(cpd, ctrl)  
}

\arguments{
  \item{cpd}{An object belonging to the class union \code{CPD}.}
  \item{ctrl}{An object of reference-class \code{Rcpp_CTRL}.}
}

\value{
  An object of reference-class \code{Rcpp_CPS}. 
}

\keyword{optimize}
