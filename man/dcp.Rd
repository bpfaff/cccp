\name{dcp}
\alias{dcp}

\title{
  Creating a member object of the reference-class \code{DCP}
}

\description{
  This function returns an object containing the definition of a convex
  program with non-linear constraints and (if provided) cone
  constraints. 
  The returned object is a member of the reference-class \code{DCP}. 
}

\usage{
dcp(x0, f0, g0, h0, cList = list(), nlfList = list(), nlgList = list(),
    nlhList = list(), A = NULL, b = NULL)
}

\arguments{
  \item{x0}{An object of class \code{vector} with dimension \eqn{n
      \times 1} for the initial values. The point \code{x0} must be in
    the domain of the nonlinear constraints.}
  \item{f0}{\code{function}: the scalar-valued convex and
    twice-differentiable objective function (its first argument must be
    \sQuote{\code{x}}).}
  \item{g0}{\code{function}: the gradient function of the objective (its
    first argument must be \sQuote{\code{x}}); returning a vector.}
  \item{h0}{\code{function}: the Hessian function of the objective (its
    first argument must be \sQuote{\code{x}}); returning a matrix.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either S4-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.} 
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
  \item{A}{An object of class \code{matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
}

\value{
  An object belonging to the reference-class \code{DCP}. 
}

\keyword{optimize}
