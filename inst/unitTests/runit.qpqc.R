##
## Unit testing of Quadratic Program with NNO and SOC constraints
test.QPIC <- function(){
    adat <- c(0.3, -0.4, -0.2, -0.4, 1.3,
              0.6, 1.2, -1.7, 0.3, -0.3,
              -0.3, 0.0, 0.6, -1.2, -2.0)
    A <- matrix(adat, nrow = 5, ncol = 3)
    b <- c(1.5, 0.0, -1.2, -0.7, 0.0)
    P <- crossprod(A)
    q <- -crossprod(A, b)
    G = -diag(3)
    h = rep(0, 3)
    nno1 <- nnoc(G = G, h = h)
    F = diag(3)
    g = rep(0, 3)
    d = rep(0, 3)
    f = 1
    soc1 <- socc(F = F, g = g, d = d, f = f)
    ans <- cccp(P = P, q = q, cList = list(nno1, soc1))
    checkTrue(ans$status == "optimal")
    return()
}
