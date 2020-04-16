\name{ctrl}
\alias{ctrl}

\title{
  Creating objects of reference-class CTRL
}

\description{
  This function creates an object of reference-class \code{CTRL} which
  contains optimization parameters, \emph{e.g.} the maximum number of
  iterations. 
}

\usage{
ctrl(maxiters = 100L, abstol = 1e-06, reltol = 1e-06,
     feastol = 1e-06, stepadj = 0.95, beta = 0.5, trace = TRUE)
}

\arguments{
  \item{maxiters}{\code{integer}, the maximum count of iterations.}
  \item{abstol}{\code{numeric}, the absolute level for convergence to be
    achieved.}
  \item{reltol}{\code{numeric}, the relative level for convergence to be
    achieved.}
  \item{feastol}{\code{numeric}, the feasable level for convergence to be
    achieved.}
  \item{stepadj}{\code{numeric}, step size adjustment in combined step.}
  \item{beta}{\code{numeric}, parameter in backtracking line search.}
  \item{trace}{\code{logical}, if \code{TRUE} (the default), the
    solver's progress during the iterations is shown.}
}

\value{
  An object of reference-class \code{CTRL}.
}

\note{
  Either \code{abstol} or \code{reltol} can be set to a negative real
  number. \code{feastol} must be greater than zero.
}

\seealso{
  \code{\linkS4class{Rcpp_CTRL}}
}

\keyword{optimize}
