##
## Function for creating 'NLFC' objects 
nlfc <- function(G, h){
    h <- new(NLFV, u = matrix(h), dims = nrow(G))
    new(NLFC, G = G, h = h, dims = nrow(G))    
}
