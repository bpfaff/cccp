\name{cccp}
\alias{cccp}

\title{
  Solving linear and quadratic programs with cone constraints
}

\description{
  This function is the main function for defining and solving convex
  problems in the form of either linear or quadratic programs with cone
  constraints. 
}

\usage{
cccp(P = NULL, q = NULL, A = NULL, b = NULL, cList = list(), optctrl = ctrl())
}

\arguments{
  \item{P}{An object of class \code{matrix} with dimension \eqn{N \times
      N} or \code{NULL}.}
  \item{q}{An object of class \code{vector} with dimension \eqn{N \times
      1} or \code{NULL}.}
  \item{A}{An object of class \code{matrix} with dimension \eqn{p \times N}.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times 1}.}
  \item{cList}{A \code{list} object containing the cone constraints.
    Elements must be of either S4-class \code{NNOC}, or \code{SOCC}, or
    \code{PSDC}.}
  \item{optctrl}{An object of S4-class \code{CTRL}.}
}

\value{
  An object of class \code{Rcpp_CPS}. 
}

\keyword{optimize}
