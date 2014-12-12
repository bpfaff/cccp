\name{psdc}
\alias{psdc}

\title{
  Definition of positive semidefinite cone inequality constraints
}

\description{
  This function is the interface to the reference class \code{PSDC} for
  creating positive semidefinite cone constraints.
}

\usage{
psdc(Flist, F0)
}

\arguments{
  \item{Flist}{Object of class \code{"list"}: A list with the matrices
  appearing on the left-hand side of the matrix inequality.} 
  \item{F0}{Object of class \code{"matrix"}: The matrix appearing on the
  righthand-side.}
}

\details{
  A psd-cone constraint is given as \eqn{\sum_{i = 1}^n x_i F_i \le
  F_0}. The matrix \eqn{G} is created as \eqn{G = [\textrm{vech}(F_1) |
  \ldots | \textrm{vech}(F_n)]} and the vector \eqn{h} is constructed as
  \eqn{h = [\textrm{vech}(F_0)]}.  
}

\value{
  List with elements: \code{conType}, \code{G} and \code{h}.
}

\keyword{classes}
