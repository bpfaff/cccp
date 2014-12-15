##
## Demo for solving an unconstrained QP
##
## Creating objects for QP
n <- 4L
M <- matrix(rnorm(n^2), nrow = n, ncol = n)
P <- crossprod(M)
q <- rnorm(n)
## Defining and solving QP
cpd <- dqp(P = P, q = q)
ctl <- ctrl()
ans <- cpd$cps(ctl)
ans
getx(ans)
