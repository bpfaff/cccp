##
## show-methods for reference objects
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
setMethod("show", "Rcpp_CTRL", function(object){
    cat("Control parameters used in optimization:\n\n")
    cat(paste("Maximum iterations:\t", object$maxiters,"\n"))
    cat(paste("Absolute tolerance:\t", object$abstol,"\n"))
    cat(paste("Relative tolerance:\t", object$reltol,"\n"))
    cat(paste("Feasible tolerance:\t", object$feastol,"\n"))
    cat(paste("Tracing progress:\t", object$trace,"\n"))
})
setMethod("show", signature = "Rcpp_DQP", function(object){
    title <- paste("* Definition of Quadratic Program *")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", ncol(object$P), "\n"))
    cat(paste("Count of equality constraints:", nrow(object$A), "\n"))
    countcc <- length(object$cList)
    cat(paste("Count of cone constraints:", countcc, "\n"))
    cc <- unlist(lapply(object$cList, function(x) class(x)))
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "Rcpp_NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "Rcpp_SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "Rcpp_PSDC")), "\n"))
    cat("\n")
})
setMethod("show", signature = "Rcpp_CPS", function(object){
    title <- "* Solution of Convex Program *"
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(title, "\n")
    cat(row, "\n")
    cat("\n")
    state <- object$state
    cat(paste("Value of primal objective:", signif(state["pobj"]), "\n"))
    if(!is.na(state["dobj"])){
        cat(paste("Value of dual objective:", signif(state["dobj"]), "\n"))
    }
    if(!is.na(state["dgap"])){
        cat(paste("Value of duality gap:", signif(state["dgap"]), "\n"))
    }
    if(!is.na(state["rdgap"])){
        cat(paste("Value of relative duality gap:", signif(state["rdgap"]), "\n"))
     }
    if(!is.na(state["certp"])){
        cat(paste("Certificate of primal infeasibility:", signif(state["certp"]), "\n"))
    }
    if(!is.na(state["certd"])){
        cat(paste("Certificate of dual infeasibility:", signif(state["certd"]), "\n"))
    }
    if(!is.na(state["pslack"])){
        cat(paste("Value of smallest primal slack:", signif(state["pslack"]), "\n"))
    }
    if(!is.na(state["dslack"])){
        cat(paste("Value of smallest dual slack:", signif(status["dslack"]), "\n"))
    }
    cat(paste("Status of solution:", object$status, "\n"))
    cat(paste("Count of iterations:", object$niter, "\n\n"))
    cat("Solutions are contained in 'pdv'.\n")
    cat("Use 'getx()', 'gety()', 'gets()' and 'getz()', respectively.\n")
})
