##
## Function for creating 'NNOC' objects 
nnoc <- function(G, h){
    h <- new(NNOV, u = matrix(h), dims = nrow(G))
    new(NNOC, G = G, h = h, dims = nrow(G))    
}
