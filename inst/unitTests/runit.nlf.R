##
## Unit testing of NLF
test.NLFV <- function(){
    n <- 10L
    unlf <- matrix(runif(n), nrow = n)
    vnlf <- matrix(runif(n), nrow = n)
    ## reference objects
    urc <- new(NLFV, u = unlf, dims = n)
    vrc <- new(NLFV, u = vnlf, dims = n)
    ## checking member functions of NLFV
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
