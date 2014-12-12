\name{nnoc}
\alias{nnoc}

\title{
  Definition of linear inequality constraints
}

\description{
  This function is the interface to the reference class \code{NNOC} for
  creating linear constraints.
}

\usage{
nnoc(G, h)
}

\arguments{
  \item{G}{Object of class \code{"matrix"}: A \eqn{(m \times n)} matrix
  containing the coefficients of the lefthand-side linear inequality
  constraints.} 
  \item{h}{Object of class \code{NNOV}: A \eqn{(m \times 1)} vector
  containing the coefficients of the righthand-side linear inequality
  constraints as slot \code{u}. 
}
}

\value{
  List with elements: \code{conType}, \code{G} and \code{h}.
}

\keyword{classes}
