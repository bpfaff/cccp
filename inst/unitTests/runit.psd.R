##
## Unit testing of PSD
test.PSDV <- function(){
    n <- 4L
    upsd <- crossprod(matrix(runif(n^2), ncol = n, nrow = n))
    dim(upsd) <- c(n^2, 1)
    vpsd <- crossprod(matrix(runif(n^2), ncol = n, nrow = n))
    dim(vpsd) <- c(n^2, 1)
    ## reference objects
    urc <- new(PSDV, u = upsd, dims = n)
    vrc <- new(PSDV, u = vpsd, dims = n)
    ## checking member functions of PSDV
    ans1 <- urc$udot(vrc)
    ans2 <- urc$uone()
    ans3 <- urc$uprd(vrc)
    ans4 <- urc$uinv(vrc)
    ans5 <- urc$umsa1(alpha = 2.0)
    ans6 <- urc$ntsc(vrc)
    ans7 <- urc$umss()
    checkTrue(ans2$dims == n)
    checkTrue(ans3$dims == n)
    checkTrue(ans4$dims == n)
    checkTrue(ans5$dims == n)
    checkTrue(nrow(ans6$r) == n)
    checkTrue(nrow(ans6$rti) == n)
    checkTrue(nrow(ans6$lambda$u) == n^2)
    checkTrue(ans6$lambda$dims == n)
    checkTrue(is.double(ans1))
    checkTrue(is.double(ans7))
    checkEqualsNumeric(sum(ans2$u), n)
    return()
}
