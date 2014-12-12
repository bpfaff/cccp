##
## Function for creating 'PSDC' objects 
psdc <- function(Flist, F0){
    m <- nrow(Flist[[1]])
    G <- do.call("cbind", lapply(Flist, function(x) c(x)))
    h <- drop(F0)
    return(list(conType = "PSDC", G = G, h = h))
}
