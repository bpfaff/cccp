\name{dlp}
\alias{dlp}

\title{
  Creating a member object of the reference-class \code{DLP}
}

\description{
  This function returns an object containing the definition of a cone
  constrained linear program. The returned object is a member of the
  reference-class \code{DLP}. 
}

\usage{
dlp(q, A = NULL, b = NULL, cList = list())
}

\arguments{
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
  An object belonging to the reference-class \code{DLP}. 
}

\keyword{optimize}
