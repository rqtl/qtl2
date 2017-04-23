# fit a single-QTL model at a single position, adjusting for polygenes with an LMM
#
fit1_pg <-
    function(genoprobs, pheno, kinship,
             addcovar=NULL, nullcovar=NULL, intcovar=NULL,
             contrasts=NULL, se=FALSE,
             hsq=NULL, reml=TRUE, tol=1e-12)
{
    stopifnot(tol > 0)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(nullcovar) && !is.matrix(nullcovar))
        nullcovar <- as.matrix(nullcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)
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

    # genoprobs is a matrix?
    if(!is.matrix(genoprobs))
        stop("genoprobs should be a matrix, individuals x genotypes")

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # check that kinship matrix is square with same IDs
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

    # multiply kinship matrix by 2; rest is using 2*kinship
    # see Almasy & Blangero (1998) https://doi.org/10.1086/301844
    kinship <- double_kinship(kinship)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, kinshipIDs,
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
    genoprobs <- genoprobs[ind2keep,,drop=FALSE]
    pheno <- pheno[ind2keep]
    if(!is.null(addcovar)) addcovar <- addcovar[ind2keep,,drop=FALSE]
    if(!is.null(nullcovar)) nullcovar <- nullcovar[ind2keep,,drop=FALSE]
    if(!is.null(intcovar)) intcovar <- intcovar[ind2keep,,drop=FALSE]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # eigen decomposition of kinship matrix
    if(!did_decomp)
        kinship <- decomp_kinship(kinship[ind2keep, ind2keep])

    # estimate hsq if necessary
    if(is.null(hsq)) {
        nullresult <- calc_hsq_clean(kinship, as.matrix(pheno), cbind(addcovar, nullcovar),
                                     NULL, FALSE, reml, cores=1, check_boundary=TRUE, tol)
        hsq <- as.numeric(nullresult$hsq)
    }

    # eigen-vectors and weights
    eigenvec <- kinship$vectors
    weights <- 1/sqrt(hsq*kinship$values + (1-hsq))

    # fit null model
    fit0 <- fit1_pg_addcovar(cbind(rep(1, length(pheno)), addcovar, nullcovar),
                             pheno,
                             matrix(ncol=0, nrow=length(pheno)),
                             eigenvec, weights, se, tol)

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs %*% contrasts

    if(is.null(intcovar)) { # just addcovar
        if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
        fitA <- fit1_pg_addcovar(genoprobs, pheno, addcovar, eigenvec, weights, se, tol)
    }
    else {                  # intcovar
        fitA <- fit1_pg_intcovar(genoprobs, pheno, addcovar, intcovar,
                                 eigenvec, weights, se, tol)
    }

    # lod score
    n <- length(pheno)
    lod <- (n/2)*log10(fit0$rss/fitA$rss)

    # residual SDs using 1/n
    sigsq0 <- fit0$sigma^2/n*fit0$df
    sigsqA <- fitA$sigma^2/n*fitA$df

    # individual contributions to the lod score
    ind_lod <- 0.5*(fit0$resid^2/sigsq0 - fitA$resid^2/sigsqA + log(sigsq0) - log(sigsqA))/log(10)
    names(ind_lod) <- names(pheno)

    # names of coefficients
    coef_names <- scan1coef_names(genoprobs, addcovar, intcovar)

    if(se) # results include standard errors
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names),
                    SE=stats::setNames(fitA$SE, coef_names)))
    else
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names)))
}
