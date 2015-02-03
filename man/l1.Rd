\name{l1}
\alias{l1}

\title{
  Minimizing L1-norm
}

\description{
  This function minimizes a L1-norm of the form \eqn{||P u - q||_1},
  whereby \eqn{P} is a \eqn{(m \times n)} matrix and \eqn{q} is a \eqn{m
  \times 1} vector. This function is wrapper function for invoking the 
  \code{cps}-method of Linear Programs.
}

\usage{
l1(P, q = NULL, optctrl = ctrl())
}

\arguments{
  \item{P}{\code{matrix} of dimension \eqn{m \times n}.}
  \item{q}{\code{vector} of length \eqn{m}.}
  \item{optctrl}{An object of S4-class \code{Rcpp_CTRL}.}
}

\value{
  An object of S4-class \code{Rcpp_CPS}.
}

\keyword{optimize}
