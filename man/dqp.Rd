\name{dqp}
\alias{dqp}

\title{
  Creating a member object of the reference-class \code{DQP}
}

\description{
  This function returns an object containing the definition of a cone
  constrained quadratic program. The returned object is a member of the
  reference-class \code{DQP}. 
}

\usage{
dqp(P, q, A = NULL, b = NULL, cList = list())
}

\arguments{
  \item{P}{An object of class \code{matrix} with dimension \eqn{n \times n}.}
  \item{q}{An object of class \code{vector} with dimension \eqn{n \times 1}.}
  \item{A}{An object of class \code{matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either reference-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.}
}

\value{
  An object belonging to the reference-class \code{DQP}. 
}

\keyword{optimize}
