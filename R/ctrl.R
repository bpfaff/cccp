##
## Function for creating 'CTRL' objects 
ctrl <- function(maxiters = 100L, abstol = 1e-7, reltol = 1e-6, feastol = 1e-7,
                 stepadj = 0.95, alpha = 0.01, beta = 0.5, maxreliter = 8L,
                 trace = TRUE){

    if(!is.integer(maxiters)){
        stop("\nThe count of maximal iterations must be an integer.\n")
    }
    if(maxiters < 1){
        stop("\nThe count of maximal iterations must be positive and greater or equal to one.\n")
    }
    if(!is.null(dim(abstol)) | length(abstol) > 1){
        stop("\nThe absolute tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(reltol)) | length(reltol) > 1){
        stop("\nThe relative tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(feastol)) | length(feastol) > 1){
        stop("\nThe feasabile tolerance for convergence must be a real scalar.\n")
    }
    if(abstol < 0 & reltol < 0){
        stop("\nAt least one of 'reltol' and 'abstol' must be positive.\n")
    }
    if(feastol <= 0){
        stop("\nThe convergence criteria for feasability must be positive.\n")
    }
    if(stepadj <= 0 || stepadj > 1.0){
        stop("\nThe step-size adjustment must be in the interval: (0, 1].\n")
    }
    if(!is.integer(maxreliter)){
        stop("\nThe count of maximum relaxed iterations must be an integer.\n")
    }
    if(maxreliter < 0){
        return("\nThe count of maximum relative iterations must be greater or equal than zero.\n")
    }
    
    new(CTRL, list(
        maxiters = maxiters,
        abstol = abstol,
        reltol = reltol,
        feastol = feastol,
        stepadj = stepadj,
        alpha = alpha,
        beta = beta,
        maxreliter = maxreliter,
        trace = as.logical(trace)[1])
        )
}
