##
## Function for creating 'NNOC' objects 
nnoc <- function(G, h){
    return(list(conType = "NNOC", G = G, h = h, dims = nrow(G)))
}
