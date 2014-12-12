\name{socc}
\alias{socc}

\title{
  Definition of second-oder cone inequality constraints
}

\description{
  This function is the interface to the reference class \code{SOCC} for
  creating second-oder cone constraints.
}

\usage{
socc(F, g, d, f)
}

\arguments{
  \item{F}{Object of class \code{"matrix"}: The matrix appearing in the
  norm-expression on the left-hand side of a second-order cone constraint.} 
  \item{g}{Object of class \code{"numeric"}: The vector appearing in the
  norm-expression on the left-hand side of a second-order cone constraint.}
  \item{d}{Object of class \code{"numeric"}: The vector appearing on the
  right-hand side of a second-order cone constraint.} 
  \item{f}{Object of class \code{"numeric"}: The scalar appearing on the
  right-hand side of a second-order cone constraint.} 
}

\details{
  A second-order cone constraint is given as \eqn{|| Fx + g ||_2 \le d'x
  + f}. The matrix \eqn{G} is created as \eqn{G = [-d, -F]} and the
  vector \eqn{h} is constructed as \eqn{h = [f, g]}. 
}

\value{
  List with elements: \code{conType}, \code{G} and \code{h}.
}

\keyword{classes}
