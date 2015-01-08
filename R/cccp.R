##
## Main function for defining and solving linear and quadratic programs with cone constraints
cccp <- function(P = NULL, q = NULL, A = NULL, b = NULL, cList = list(),
                 optctrl = ctrl()){

    if(is.null(P) & is.null(q)){
        stop("At least, 'P' or 'q' must be provided for quadratic or linear objective.\n")
    }
    if(is.null(P)){
        stop("Only Quadratic problems implemented so far.\n")
    }
    cpd <- dqp(P, q, A, b, cList)
    cpd$cps(optctrl)
}
