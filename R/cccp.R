##
## Main function for defining and solving linear and quadratic programs with cone constraints
cccp <- function(P = NULL, q = NULL, A = NULL, b = NULL, cList = list(),
                 x0 = NULL, nlfList = list(), nlgList = list(), nlhList = list(),
                 optctrl = ctrl()){

    if(is.null(P) & is.null(q)){
        stop("At least, 'P' or 'q' must be provided for quadratic or linear objective.\n")
    }
    if(is.null(P)){
        if(is.null(x0)){
            cpd <- dlp(q, A, b, cList)
        } else {
            cpd <- dnl(q, A, b, cList, x0, nlfList, nlgList, nlhList)
        }
    } else {
        cpd <- dqp(P, q, A, b, cList)
    }
    cps(cpd, optctrl)
}
