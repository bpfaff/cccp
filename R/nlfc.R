##
## Function for creating 'NLFC' objects 
nlfc <- function(G, h){
    return(list(conType = "NLFC", G = G, h = h, dims = nrow(G)))
}
