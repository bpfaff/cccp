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
cccp(P = NULL, q = NULL, A = NULL, b = NULL, cList = list(),
     x0 = NULL, f0 = NULL, g0 = NULL, h0 = NULL,
     nlfList = list(), nlgList = list(), nlhList = list(),
     optctrl = ctrl())
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
  \item{x0}{An object of class \code{vector} with dimension \eqn{n
      \times 1} for the initial values. The point \code{x0} must be in
    the domain of the nonlinear constraints.}
  \item{f0}{\code{function}: the scalar-valued convex and
    twice-differentiable objective function (its first argument must be
    \sQuote{\code{x}}).}
  \item{g0}{\code{function}: the gradient function of the objective (its
    first argument must be \sQuote{\code{x}}).}
  \item{h0}{\code{function}: the Hessian function of the objective (its
    first argument must be \sQuote{\code{x}}).}
  \item{nlfList}{A \code{list} object containing the nonlinear
    constraints as its elements. The functions have to be specified with
    \code{x} as their first argument and must be casted in implicit
    form, \emph{i.e.} \eqn{f(x) \le 0}.}
  \item{nlgList}{A \code{list} object containing the gradient functions
    as its elements. The functions have to be specified with \code{x} as
    their first argument.}
  \item{nlhList}{A \code{list} object containing the Hessian functions
    as its elements. The functions have to be specified with \code{x} as
    their first argument.}
  \item{optctrl}{An object of S4-class \code{Rcpp_CTRL}.}
}

\value{
  An object of class \code{Rcpp_CPS}. 
}

\keyword{optimize}
