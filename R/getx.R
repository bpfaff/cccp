##
## Value of variable 'x'
setMethod("getx", signature = "Rcpp_PDV", function(object){
    object$x
})
setMethod("getx", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    getx(pdv)
})
