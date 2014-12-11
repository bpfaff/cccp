##
## Value of variable 'y'
setMethod("gety", signature = "Rcpp_PDV", function(object){
    object$y
})
setMethod("gety", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    gety(pdv)
})
