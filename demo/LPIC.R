##
## Demo for solving a Linear Program with (linear) inequality constraints
## (Example taken from cvxopt's userguide)
##
## Creating LP
q <- c(-4, -5)
G <- matrix(c(2, 1, -1, 0,
              1, 2, 0, -1),
            nrow = 4, ncol = 2)
h <- c(3, 3, 0, 0)
nno1 <- nnoc(G = G, h = h)
## Using main function of package
ans <- cccp(q = q, cList = list(nno1), optctrl = ctrl())
ans
getx(ans)
