# Calculate QTL effects in scan along one chromosome, adjusting for polygenes with an LMM
#
scan1coef_pg <-
    function(genoprobs, pheno, kinship,
             addcovar=NULL, nullcovar=NULL, intcovar=NULL,
             weights=NULL, contrasts=NULL, se=FALSE,
             hsq=NULL, reml=TRUE, ...)
{
    # deal with the dot args
    dotargs <- list("...")
    tol <- grab_dots(dotargs, "tol", 1e-12)
    if(!is_pos_number(tol)) stop("tol should be a single positive number")
    check_extra_dots(dotargs, "tol")

    # check that the objects have rownames
    check4names(pheno, addcovar, NULL, intcovar, nullcovar)

    # force things to be matrices
    if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.numeric(addcovar)) stop("addcovar is not numeric")
    }
    if(!is.null(nullcovar)) {
        if(!is.matrix(nullcovar)) nullcovar <- as.matrix(nullcovar)
        if(!is.numeric(nullcovar)) stop("nullcovar is not numeric")
    }
    if(!is.null(intcovar)) {
        if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)
        if(!is.numeric(intcovar)) stop("intcovar is not numeric")
    }
    if(!is.null(contrasts)) {
        if(!is.matrix(contrasts)) contrasts <- as.matrix(contrasts)
        if(!is.numeric(contrasts)) stop("contrasts is not numeric")
    }

    # make sure pheno is a vector
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno) > 1)
            warning("Considering only the first phenotype.")
        rn <- rownames(pheno)
        pheno <- pheno[,1]
        names(pheno) <- rn
        if(!is.numeric(pheno)) stop("pheno is not numeric")
    }

    # genoprobs has more than one chromosome?
    if(length(genoprobs) > 1)
        warning("Using only the first chromosome, ", names(genoprobs)[1])
    genoprobs <- genoprobs[[1]]

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # make sure kinship is for a single chromosome and get IDs
    did_decomp <- is_kinship_decomposed(kinship)
    kinship <- check_kinship_onechr(kinship)
    kinshipIDs <- check_kinship(kinship, 1)

    # multiply kinship matrix by 2; rest is using 2*kinship
    # see Almasy & Blangero (1998) https://doi.org/10.1086/301844
    kinship <- double_kinship(kinship)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, kinshipIDs, weights,
                               addcovar, nullcovar, intcovar, complete.cases=TRUE)

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
    if(!is.null(nullcovar)) nullcovar <- nullcovar[ind2keep,,drop=FALSE]
    if(!is.null(intcovar)) intcovar <- intcovar[ind2keep,,drop=FALSE]
    if(!is.null(weights)) weights <- weights[ind2keep]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # square-root of weights (only if model="normal")
    weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)

    # multiply stuff by weights
    kinship <- weight_kinship(kinship, weights)
    addcovar <- weight_matrix(addcovar, weights)
    intcovar <- weight_matrix(intcovar, weights)
    nullcovar <- weight_matrix(nullcovar, weights)
    pheno <- weight_matrix(pheno, weights)

    # eigen decomposition of kinship matrix
    if(!did_decomp)
        kinship <- decomp_kinship(kinship[ind2keep, ind2keep])

    # estimate hsq if necessary
    if(is.null(hsq)) {
        nullresult <- calc_hsq_clean(Ke=kinship, pheno=as.matrix(pheno),
                                     addcovar=cbind(addcovar, nullcovar), Xcovar=NULL,
                                     weights=weights, is_x_chr=FALSE, reml=reml,
                                     cores=1, check_boundary=TRUE, tol=tol)
        hsq <- as.numeric(nullresult$hsq)
    }

    # eigen-vectors and weights
    eigenvec <- kinship$vectors
    wts <- 1/sqrt(hsq*kinship$values + (1-hsq))

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    # multiply genoprobs by weights
    genoprobs <- weight_array(genoprobs[ind2keep,,,drop=FALSE], weights)

    if(se) { # also calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            result <- scancoefSE_pg_addcovar(genoprobs, pheno, addcovar, eigenvec, wts, tol)
        }
        else {                  # intcovar
            result <- scancoefSE_pg_intcovar(genoprobs, pheno, addcovar, intcovar,
                                              eigenvec, wts, tol)
        }

        SE <- t(result$SE) # transpose to positions x coefficients
        result <- result$coef
    } else { # don't calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            result <- scancoef_pg_addcovar(genoprobs, pheno, addcovar, eigenvec, wts, tol)
        }
        else {                  # intcovar
            result <- scancoef_pg_intcovar(genoprobs, pheno, addcovar, intcovar,
                                            eigenvec, wts, tol)
        }
        SE <- NULL
    }

    result <- t(result) # transpose to positions x coefficients

    # add names
    dimnames(result) <- list(dimnames(genoprobs)[[3]],
                             scan1coef_names(genoprobs, addcovar, intcovar))
    if(se) dimnames(SE) <- dimnames(result)

    # add some attributes with details on analysis
    attr(result, "sample_size") <- length(ind2keep)
    attr(result, "SE") <- SE # include only if not NULL

    class(result) <- c("scan1coef", "scan1", "matrix")
    result
}
