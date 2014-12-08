library(cccp3)
library(RUnit)
##
## Unit testing of SOC
test.SOCV <- function(){
    n <- 10L
    usoc <- matrix(c(n + 1, runif(n - 1)), nrow = n)
    vsoc <- matrix(c(n + 1, runif(n - 1)), nrow = n)
    ## reference objects
    urc <- new(SOCV, u = usoc, dims = n)
    vrc <- new(SOCV, u = vsoc, dims = n)
    ## checking member functions of SOCV
    ans1 <- urc$udot(vrc)
    ans2 <- urc$uone()
    ans3 <- urc$uprd(vrc)
    ans4 <- urc$uinv(vrc)
    ans5 <- urc$umsa(alpha = 2.0, init = TRUE)
    ans6 <- urc$umsa(alpha = 2.0, init = FALSE)
    ans7 <- urc$ntsc(vrc)
    ans8 <- urc$umss()
    ans9 <- urc$jdot(vrc)
    checkTrue(ans2$dims == n)
    checkTrue(ans3$dims == n)
    checkTrue(ans4$dims == n)
    checkTrue(ans5$dims == n)
    checkTrue(ans6$dims == n)
    checkTrue(nrow(ans7$v) == n)
    checkTrue(is.double(ans7$beta))
    checkTrue(nrow(ans7$lambda$u) == n)
    checkTrue(ans7$lambda$dims == n)
    checkTrue(is.double(ans1))
    checkTrue(is.double(ans8))
    checkTrue(is.double(ans9))
    checkEqualsNumeric(sum(ans2$u), 1.0)
    checkEqualsNumeric(urc$u, vrc$uprd(ans4)$u)
    return()
}
