\name{nlfc}
\alias{nlfc}

\title{
  Definition of nonlinear inequality constraints
}

\description{
  This function is the interface to the reference class \code{NLFC} for
  creating nonlinear constraints.
}

\usage{
nlfc(G, h)
}

\arguments{
  \item{G}{Object of class \code{"matrix"}: A \eqn{(m \times n)} matrix
  containing the coefficients of the lefthand-side linear inequality
  constraints.} 
  \item{h}{Object of class \code{NLFV}: A \eqn{(m \times 1)} vector
  containing the coefficients of the righthand-side linear inequality
  constraints as slot \code{u}. 
}
}

\value{
  List with elements: \code{conType}, \code{G} and \code{h}.
}

\keyword{classes}
