\name{CPG}
\alias{CONEC}
\alias{CTRL}
\alias{CPG}
\alias{PDV}
\alias{DCP}
\alias{DLP}
\alias{DNL}
\alias{DQP}
\alias{CPS}
\alias{rpp}
\alias{gpp}

\title{
	Rcpp module: CPG
}

\description{
  Module for defining and solving convex programs.
}

\details{
  The module contains the following items:
  classes: \describe{
    \item{CONEC}{Class for inequality (cone) constraints.}
    \item{CTRL}{Class for control parameters used in optimizations.}
    \item{PDV}{Class for primal/dual variables.}
    \item{DCP}{Class for definition of convex programs.}
    \item{DLP}{Class for definition of linear programs.}
    \item{DNL}{Class for definition of linear programs with non-linear constraints.}
    \item{DQP}{Class for definition of quadratic programs.}
    \item{CPS}{Class for solution of convex programs.}
  }
  functions: \describe{
    \item{rpp}{Function for solving risk parity portfolios.}
    \item{gpp}{Function for solving a geometric program.}
  }
}
    
\keyword{datasets}
