# Calculate BLUPs of QTL effects in scan along one chromosome, with residual polygenic effect
scan1blup_pg <-
    function(genoprobs, pheno, kinship, addcovar=NULL,
             contrasts=NULL, se=FALSE, reml=TRUE, preserve_intercept=FALSE,
             tol=1e-12, cores=1, quiet=TRUE)
{
    stopifnot(tol > 0)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(contrasts) && !is.matrix(contrasts))
        contrasts <- as.matrix(contrasts)

    # make sure pheno is a vector
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno) > 1)
            warning("Considering only the first phenotype.")
        rn <- rownames(pheno)
        pheno <- pheno[,1]
        names(pheno) <- rn
    }

    # genoprobs has more than one chromosome?
    if(length(genoprobs$probs) > 1)
        warning("Using only the first chromosome, ", names(genoprobs)[1])
    map <- genoprobs$map[[1]]
    genoprobs <- genoprobs$probs[[1]]

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # check that kinship matrices are square with same IDs
    if(!is.null(attr(kinship, "eigen_decomp"))) { # already did decomposition
        kinshipIDs <- rownames(kinship$vectors)
        did_decomp <- TRUE
    } else {
        if(is.list(kinship)) { # if a list of length one, take the first part
            if(length(kinship) > 1) # if a list of length >1, give error
                stop("kinship should be a single matrix")
            kinship <- kinship[[1]]
        }
        kinshipIDs <- check_kinship(kinship, 1)
        did_decomp <- FALSE
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, kinshipIDs, addcovar, complete.cases=TRUE)
    if(length(ind2keep)<=2) {
        if(length(ind2keep)==0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse=":"))
    }

    if(did_decomp) { # if did decomposition already, make sure it was with exactly
        if(length(kinshipIDs) != length(ind2keep) ||
           any(sort(kinshipIDs) != sort(ind2keep)))
            stop("Decomposed kinship matrix was with different individuals")
        else
            ind2keep <- kinshipIDs # force them in same order
    }

    # omit individuals not in common
    genoprobs <- genoprobs[ind2keep,,,drop=FALSE]
    pheno <- pheno[ind2keep]
    if(!is.null(addcovar)) addcovar <- addcovar[ind2keep,,drop=FALSE]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)
    # add the intercept
    addcovar <- cbind(intercept=rep(1, length(ind2keep)),
                      addcovar)
    rownames(addcovar) <- ind2keep

    # eigen decomposition of kinship matrix
    if(!did_decomp)
        kinship <- decomp_kinship(kinship[ind2keep, ind2keep])


    # estimate hsq if necessary
    nullresult <- calc_hsq_clean(kinship, as.matrix(pheno), addcovar, NULL, FALSE,
                                 reml, cores=1, check_boundary=TRUE, tol)
    hsq <- nullresult$hsq

    # eigen-vectors and weights
    eigenvec <- kinship$vectors
    weights <- 1/sqrt(hsq*kinship$values + (1-hsq))

    # rotate genoprobs, pheno, and addcovar
    genoprobs <- matrix_x_3darray(eigenvec, genoprobs)
    addcovar <- eigenvec %*% addcovar
    pheno <- eigenvec %*% pheno

    #  multiply by weights
    genoprobs <- weighted_3darray(genoprobs, weights)
    pheno <- pheno * weights
    addcovar <- addcovar * weights

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # determine batches
    n_pos <- dim(genoprobs)[[3]]
    batches <- batch_vec(1:n_pos, ceiling(n_pos/n_cores(cores)))

    by_group_func <- function(i)
        scanblup(genoprobs[,,batches[[i]],drop=FALSE], pheno, addcovar, se, reml, preserve_intercept, tol)

    # scan to get BLUPs and coefficient estimates
    if(n_cores(cores)==1) {
        result <- scanblup(genoprobs, pheno, addcovar, se, reml, preserve_intercept, tol)
        if(se) SE <- t(result$SE)
        else SE <- NULL
        coef <- t(result$coef)
    } else {
        result <- cluster_lapply(cores, seq(along=batches), by_group_func)
        SE <- coef <- matrix(nrow=n_pos,ncol=nrow(result[[1]]$coef))
        for(i in seq(along=result)) {
            coef[batches[[i]],] <- t(result[[i]]$coef)
            if(se) SE[batches[[i]],] <- t(result[[i]]$SE)
            else SE <- NULL
        }
    }

    # add names
    if(preserve_intercept)
        coefnames <- scan1coef_names(genoprobs, addcovar, NULL)
    else
        coefnames <- scan1coef_names(genoprobs, addcovar[,-1,drop=FALSE], NULL)
    dimnames(coef) <- list(dimnames(genoprobs)[[3]], coefnames)
    if(se) dimnames(SE) <- dimnames(coef)

    result <- list(coef = coef,
                   map = map,
                   sample_size = length(ind2keep),
                   addcovar = colnames4attr(addcovar),
                   contrasts = contrasts)
    result$SE <- SE # include only if not NULL

    class(result) <- c("scan1coef", "scan1", "matrix")
    result
}
