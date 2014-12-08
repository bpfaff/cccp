##
## Unit testing of NNO
test.NNOV <- function(){
    n <- 10L
    unno <- matrix(runif(n), nrow = n)
    vnno <- matrix(runif(n), nrow = n)
    ## reference objects
    urc <- new(NNOV, u = unno, dims = n)
    vrc <- new(NNOV, u = vnno, dims = n)
    ## checking member functions of NNOV
    ans1 <- urc$udot(vrc)
    ans2 <- urc$uone()
    ans3 <- urc$uprd(vrc)
    ans4 <- urc$uinv(vrc)
    ans5 <- urc$umsa(alpha = 2.0, init = TRUE)
    ans6 <- urc$umsa(alpha = 2.0, init = FALSE)
    ans7 <- urc$ntsc(vrc)
    ans8 <- urc$umss()
    checkTrue(ans2$dims == n)
    checkTrue(ans3$dims == n)
    checkTrue(ans4$dims == n)
    checkTrue(ans5$dims == n)
    checkTrue(ans6$dims == n)
    checkTrue(nrow(ans7$dnl) == n)
    checkTrue(nrow(ans7$dnli) == n)
    checkTrue(nrow(ans7$lambda$u) == n)
    checkTrue(ans7$lambda$dims == n)
    checkTrue(is.double(ans1))
    checkTrue(is.double(ans8))
    checkEqualsNumeric(sum(ans2$u), n)
    checkEqualsNumeric(urc$u, vrc$uprd(ans4)$u)
    return()
}
