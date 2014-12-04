loadModule("NLF", TRUE)
loadModule("NNO", TRUE)
loadModule("SOC", TRUE)
loadModule("PSD", TRUE)

evalqOnLoad({
    ## show methods
    setMethod("show", "Rcpp_NLFV", function(object){
        cat("Object of reference-class: NLFV\n")
        cat(paste("Dimension:", object$dims,"\n"))
        cat("\n")
        cat("Showing head of NLF-variable\n")
        print(head(object$u))
    })
    setMethod("show", "Rcpp_NNOV", function(object){
        cat("Object of reference-class: NNOV\n")
        cat(paste("Dimension:", object$dims,"\n"))
        cat("\n")
        cat("Showing head of NNO-variable\n")
        print(head(object$u))
    })
    setMethod("show", "Rcpp_SOCV", function(object){
        cat("Object of reference-class: SOCV\n")
        cat(paste("Dimension:", object$dims,"\n"))
        cat("\n")
        cat("Showing head of SOC-variable\n")
        print(head(object$u))
    })
    setMethod("show", "Rcpp_PSDV", function(object){
        cat("Object of reference-class: PSDV\n")
        cat(paste("Dimension:", object$dims,"\n"))
        cat("\n")
        cat("Showing head of PSD-variable\n")
        print(head(object$u))
    })
})
