\name{dnl}
\alias{dnl}

\title{
  Creating a member object of the reference-class \code{DNL}
}

\description{
  This function returns an object containing the definition of a linear
  program with non-linear constraints and (if provided) cone
  constraints. 
  The returned object is a member of the reference-class \code{DNL}. 
}

\usage{
dnl(q, A = NULL, b = NULL, cList = list(),
    x0, nlfList = list(), nlgList = list(), nlhList = list())
}

\arguments{
  \item{q}{\code{vector} of length \eqn{n} for the coefficients in the
    objective.}
  \item{A}{An object of class \code{matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either S4-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.}
  \item{x0}{An object of class \code{vector} with dimension \eqn{n
      \times 1} for the initial values. The point \code{x0} must be in
    the domain of the nonlinear constraints.}
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
}

\value{
  An object belonging to the reference-class \code{DNL}. 
}

\keyword{optimize}
