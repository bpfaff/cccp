\name{rp}
\alias{rp}

\title{
  Risk-parity optimization
}

\description{
  This function determines a risk-parity solution of a long-only
  portfolio with a budget-constraint.
}

\usage{
rp(x0, P, mrc, optctrl = ctrl())
}

\arguments{
  \item{x0}{\code{matrix} of dimension \eqn{n \times 1}; starting values.}
  \item{P}{\code{matrix} of dimension \eqn{n \times n}; dispersion matrix.}
  \item{mrc}{\code{matrix} of dimension \eqn{n \times 1}; the marginal
    risk contributions.}
  \item{optctrl}{An object of S4-class \code{Rcpp_CTRL}.}
}

\value{
  An object of S4-class \code{Rcpp_CPS}.
}

\references{
  Spinu, F. (2013), An Algorithm for Computing Risk Parity Weights,
  SSRN, \emph{OMERS Capital Markets}, July 2013.
}

\keyword{optimize}
