##
## Function for creating an object of reference-class 'DQP'
dqp <- function(P, q, A = NULL, b = NULL, cList = list()){
   n <- ncol(P)
   if(is.null(A)){
       A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    }
    if(is.null(b)){
        b <- numeric()
    } 
    if(length(cList) > 0){
        coneclasses <- unlist(lapply(cList, class))
        if(!all(coneclasses %in% c("Rcpp_NNOC", "Rcpp_SOCC", "Rcpp_PSDC"))){
            stop("List elements of cone constraints must be of either class:\n'Rcpp_NNOC', or 'Rcpp_SOCC', or 'Rcpp_PSDC'.\n")
        }
    }
   ans <- new(DQP,
              P = P,
              q = q,
              A = A,
              b = b,
              cList = cList)
   return(ans)
}
