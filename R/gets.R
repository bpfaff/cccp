##
## Value of variable 's'
setMethod("gets", signature = "Rcpp_PDV", function(object){
    lapply(object$s, function(s) s$u)
})
setMethod("gets", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    gets(pdv)
})
