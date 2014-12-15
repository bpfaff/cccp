##
## Demo for solving an unconstrained QP
##
## Creating QP
## Objective
n <- 4L
M <- matrix(rnorm(n^2), nrow = n, ncol = n)
P <- crossprod(M)
q <- rnorm(n)
## Using main function of package
cpd <- dqp(P = P, q = q)
ctl <- ctrl()
cpd$cps(ctl)
