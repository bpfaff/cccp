##
## Function for creating 'PSDC' objects 
psdc <- function(Flist, F0){
    m <- nrow(Flist[[1]])
    G <- do.call("cbind", lapply(Flist, function(x) c(x)))
    h <- new(PSDV, u = matrix(as(F0, "vector"), ncol = 1), dims = m)
    new(PSDC, G = G, h = h, dims = m)        
}
