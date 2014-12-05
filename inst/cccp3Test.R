##
## Loading development version of package cccp
## (module-based implementation with C++ classes)
## 
library(cccp3)
##
## Creating NLFV objects
##
n <- 10L
unlf <- matrix(runif(n), nrow = n)
vnlf <- matrix(runif(n), nrow = n)
## reference objects
urc <- new(NLFV, u = unlf, dims = n)
vrc <- new(NLFV, u = vnlf, dims = n)
## udot-methods
(rcobj <- urc$udot(vrc))
## uone-methods
(rcobj <- urc$uone())
## uprd-methods
(rcobj <- urc$uprd(vrc))
## uinv-methods
(rcobj <- urc$uinv(vrc))
## umsa-methods
(rcobj <- urc$umsa(alpha = 2.0, init = TRUE))
(rcobj <- urc$umsa(alpha = 2.0, init = FALSE))
## ntsc-method
(rcobj <- urc$ntsc(vrc))
rcobj$dnl
rcobj$dnli
rcobj$lambda
##
## Creating NNOV objects
##
n <- 10L
unno <- matrix(runif(n), nrow = n)
vnno <- matrix(runif(n), nrow = n)
## reference objects
urc <- new(NNOV, u = unno, dims = n)
vrc <- new(NNOV, u = vnno, dims = n)
## udot-methods
(rcobj <- urc$udot(vrc))
## uone-methods
(rcobj <- urc$uone())
## uprd-methods
(rcobj <- urc$uprd(vrc))
## uinv-methods
(rcobj <- urc$uinv(vrc))
## umsa-methods
(rcobj <- urc$umsa(alpha = 2.0, init = TRUE))
## ntsc-method
(rcobj <- urc$ntsc(vrc))
rcobj$d
rcobj$di
rcobj$lambda
##
## Creating SOCV objects
##
n <- 10L
usoc <- matrix(c(n + 1, runif(n - 1)), nrow = n)
vsoc <- matrix(c(n + 1, runif(n - 1)), nrow = n)
## reference objects
urc <- new(SOCV, u = usoc, dims = n)
vrc <- new(SOCV, u = vsoc, dims = n)
## udot-methods
(rcobj <- urc$udot(vrc))
## uone-methods
(rcobj <- urc$uone())
## uprd-methods
(rcobj <- urc$uprd(vrc))
## uinv-methods
(rcobj <- urc$uinv(vrc))
## umsa-methods
(rcobj <- urc$umsa(alpha = 2.0, init = TRUE))
## ntsc-method
(rcobj <- urc$ntsc(vrc))
rcobj$v
rcobj$beta
rcobj$lambda
##
## Creating PSDV objects
##
n <- 10L
upsd <- crossprod(matrix(rnorm(n^2), ncol = n, nrow = n))
dim(upsd) <- c(n^2, 1)
vpsd <- crossprod(matrix(rnorm(n^2), ncol = n, nrow = n))
dim(vpsd) <- c(n^2, 1)
## reference objects
urc <- new(PSDV, u = upsd, dims = n)
vrc <- new(PSDV, u = vpsd, dims = n)
## udot-methods
(rcobj <- urc$udot(vrc))
## uone-methods
(rcobj <- urc$uone())
## uprd-methods
(rcobj <- urc$uprd(vrc))
## uinv-methods
(rcobj <- urc$uinv(vrc))
## umsa-methods
(rcobj <- urc$umsa(alpha = 2.0, init = TRUE))
## ntsc-method
