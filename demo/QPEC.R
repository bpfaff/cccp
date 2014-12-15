##
## Demo for solving an unconstrained QP
##
## Creating objects for QP
P <- 2 * matrix(c(2, .5, .5, 1), nrow = 2, ncol = 2)
q <- c(1.0, 1.0)
A <- matrix(c(1.0, 1.0), nrow = 1, ncol = 2)
b <- 1.0
## Defining and solving QP
cpd <- dqp(P = P, q = q, A = A, b = b)
ctl <- ctrl()
ans <- cpd$cps(ctl)
ans
getx(ans)
