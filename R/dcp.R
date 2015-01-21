##
## Function for creating an object of reference-class 'DCP'
dcp <- function(x0, f0, g0, h0, cList = list(), nlfList = list(), nlgList = list(), nlhList = list(), A = NULL, b = NULL){
    x0 <- as.matrix(x0)
    n <- nrow(x0)
    K <- length(cList)
    mnl <- length(nlfList)
    ## Checking whether x0 is in the domain of nonlinear objective
    f0Dom <- is.nan(f0(x0))
    if(f0Dom){
        stop("Initial point 'x0' is not in the domain of nonlinear objective 'f0'.\n")
    }
    ## Creating list-object of non-linear objective, its Gradient and Hessian functions
    oList <- list(f0, g0, h0)
    ##
    ## Checking provided non-linear constraints (if applicable)
    ##
    if(mnl > 0){
        if(!all(unlist(lapply(nlfList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlfList' are functions.\n")
        }
        if(!all(unlist(lapply(nlgList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlgList' are functions.\n")
        }
        if(!all(unlist(lapply(nlhList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlhList' are functions.\n")
        }
        fDom <- unlist(lapply(nlfList, function(fcc) fcc(x0)))
        idxnan <- which(is.nan(fDom))
        if(any(idxnan)){
            stop(paste("Initial point 'x0' is not in the domain of nonlinear convex constraint(s): ", idxnan, ".\n", sep = ""))
        }
        if(length(nlfList) != length(nlgList)){
            stop("Length of lists for nonlinear functions and gradient functions do differ.\n")
        }
        if(length(nlfList) != length(nlhList)){
            stop("Length of lists for nonlinear functions and Hessian functions do differ.\n")
        }
        ## Creating list-object of non-linear constraints, their Gradient and Hessian functions
        nList <- list(nlfList, nlgList, nlhList)
        ## Creating objects related to NLFC
        Gnl <- matrix(0, nrow = mnl, ncol = n)
        hnl <- matrix(0, nrow = mnl, ncol = 1)
    } else {
        nList <- list()
    }
    ##
    ## Checking cone constraints
    ##
    if(K > 0){
        cone <- unlist(lapply(cList, function(x) x[["conType"]]))
        if(!all(cone %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be either created by calls to:\n'nnoc()', or 'socc()', or 'psdc()'.\n")
        }
        GList <- lapply(cList, function(x) x[["G"]])
        hList <- lapply(cList, function(x) x[["h"]])
        if(mnl > 0){
            cone <- c("NLFC", cone)
            GList <- c(list(Gnl), GList)
            hList <- c(list(hnl), hList)
            dims <- c(mnl, as.integer(unlist(lapply(cList, function(x) x[["dims"]]))))
            K <- K + 1L
        }
        G <- do.call("rbind", GList)
        h <- do.call("rbind", hList)
        ridx <- cumsum(unlist(lapply(GList, nrow)))
        sidx <- cbind(c(0, ridx[-length(ridx)]), ridx - 1)
        dims <- as.integer(unlist(lapply(cList, function(x) x[["dims"]])))
        cList <- new(CONEC, cone, G, h, sidx, dims, K, n)
    } else if(mnl > 0){ ## case: no cone constraints, but nonlinear constraints
        cList <- new(CONEC, "NLFC", Gnl, hnl, mnl, 1L, n)
    } else { ## case: neither cone nor nonlinear constraints
        cList <- new(CONEC, as.integer(n))
    }
    ##
    ## Checking equality constraints
    ##
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    } 
    if(is.null(b)){
        b <- matrix(0, nrow = 0, ncol = 1)
    }
    if(is.null(dim(b))){
        b <- matrix(b, ncol = 1)
    } 
    ## checking whether x0 satisfies equality constraints
    if(!is.null(A)){
        eq <- identical(as.vector(A %*% x0 - b), rep(0, nrow(A)))
        if(!eq){
            stop("Initial point 'x0' does not satisfy equality constraints.\n")
        }
    }
    
    ans <- new(DCP,
               x0 = as.matrix(x0),
               cList = cList,
               oList = oList,
               nList = nList,
               A = A,
               b = b
               )
    return(ans)
}
