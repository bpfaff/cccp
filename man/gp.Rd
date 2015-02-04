\name{gp}
\alias{gp}

\title{
  Geometric program
}

\description{
  This function solves a geometric program.
}

\usage{
gp(F0, g0, FList = list(), gList = list(), nno = NULL,
   A = NULL, b = NULL, optctrl = ctrl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F0}{Matrix in the objective function.}
  \item{g0}{Matrix in the objective function (affine terms).}
  \item{FList}{List of matrices in posinomial functions.}
  \item{gList}{List of matrices in posinomial functions (affine terms).}
  \item{nno}{Object created by a call to \code{nnoc()}.}
  \item{A}{Lefthand-side matrix of equality constraints.}
  \item{b}{Lefthand-side matrix of equality constraints.}
  \item{optctrl}{Object of reference class \sQuote{Rcpp_CTRL}, created
    by a call to \code{ctrl()}.}
}

\details{
  Solves a geometric program casted in its epigraph form.
}

\value{
  An object of S4-class \code{Rcpp_CPS}.
}

\references{
  Boyd, S., Kim, S.-J., Vandenberghe, L. and A. Hassibi (2007), A tutorial on
  geometric programming, \emph{Optim Eng}, Educational Section,
  \bold{8}:67--127, Springer. 
}

\keyword{optimize}
