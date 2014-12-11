##
## Value of variable 'z'
setMethod("getz", signature = "Rcpp_PDV", function(object){
    lapply(object$z, function(s) z$u)
})
setMethod("getz", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    getz(pdv)
})
