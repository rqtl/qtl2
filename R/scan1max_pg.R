# Genome scan with a single-QTL and linear mixed model
#
# called by scan1max() when kinship() is provided.
scan1max_pg <-
    function(genoprobs, pheno, kinship, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, hsq=NULL,
             by_chr=FALSE, cores=1, ...)
{
    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    if(!is_pos_number(tol)) stop("tol should be a single positive number")
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    if(!is.null(max_batch) && !is_pos_number(max_batch)) stop("max_batch should be a single positive integer")
    check_boundary <- grab_dots(dotargs, "check_boundary", TRUE)
    check_extra_dots(dotargs, c("tol", "intcovar_method", "check_boundary", "quiet", "max_batch"))

    # check that the objects have rownames
    check4names(pheno, addcovar, Xcovar, intcovar)

    # force things to be matrices
    if(!is.matrix(pheno)) {
        pheno <- as.matrix(pheno)
        if(!is.numeric(pheno)) stop("pheno is not numeric")
    }
    if(is.null(colnames(pheno))) # force column names
        colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
    if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.numeric(addcovar)) stop("addcovar is not numeric")
    }
    if(!is.null(Xcovar)) {
        if(!is.matrix(Xcovar)) Xcovar <- as.matrix(Xcovar)
        if(!is.numeric(Xcovar)) stop("Xcovar is not numeric")
    }
    if(!is.null(intcovar)) {
        if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)
        if(!is.numeric(intcovar)) stop("intcovar is not numeric")
    }

    # check that kinship matrices are square with same IDs
    kinshipIDs <- check_kinship(kinship, length(genoprobs))

    # multiply kinship matrix by 2; rest is using 2*kinship
    # see Almasy & Blangero (1998) https://doi.org/10.1086/301844
    kinship <- double_kinship(kinship)

    # take square-root of weights
    weights <- sqrt_weights(weights)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
                               kinshipIDs, weights, complete.cases=TRUE)
    ind2keep <- get_common_ids(ind2keep, pheno[rowSums(is.finite(pheno)) > 0,,drop=FALSE])
    if(length(ind2keep)<=2) {
        if(length(ind2keep)==0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse=":"))
    }

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # drop things from Xcovar that are already in addcovar
    Xcovar <- drop_xcovar(addcovar, Xcovar, tol)

    # batch phenotypes by missing values
    phe_batches <- batch_cols(pheno[ind2keep,,drop=FALSE], max_batch)

    # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
    genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    n_chr <- length(genoprobs)

    # to contain the results
    result <- matrix(nrow=n_chr, ncol=ncol(pheno))
    dimnames(result) <- list(names(genoprobs), colnames(pheno))

    # number of chr to consider under null
    if(is_kinship_list(kinship)) n_null_chr <- length(kinship)
    else if(!is.null(Xcovar)) n_null_chr <- any(is_x_chr) + any(!is_x_chr)
    else n_null_chr <- 1

    if(is.null(hsq)) {
        hsq <- matrix(nrow=n_null_chr, ncol=ncol(pheno))
        dimnames(hsq) <- hsq_dimnames(kinship, Xcovar, is_x_chr, pheno)
        estimate_hsq <- TRUE
    } else {
        if(length(hsq) != n_null_chr * ncol(pheno)) {
            stop("hsq should be NULL or a matrix of size ", n_null_chr, " x ", ncol(pheno))
        }
        if(!is.matrix(hsq)) hsq <- as.matrix(hsq, ncol=ncol(pheno))
        if(nrow(hsq) != n_null_chr) {
            stop("hsq should be NULL or a matrix of size ", n_null_chr, " x ", ncol(pheno))
        }
        if(any(hsq < 0 | hsq > 1)) {
            stop("hsq values should be between 0 and 1.")
        }
        estimate_hsq <- FALSE
    }
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)

    chr_index <- rep(seq_len(n_chr), sapply(genoprobs, function(a) dim(a)[3]))

    # loop over batches of phenotypes with the same pattern of NAs
    for(batch in seq_along(phe_batches)) {

        # info about batch
        omit <- phe_batches[[batch]]$omit # ind to omit
        phecol <- phe_batches[[batch]]$cols # phenotype columns in batch

        # individuals to keep in this batch
        these2keep <- ind2keep
        if(length(omit)>0) these2keep <- ind2keep[-omit]
        n[phecol] <- length(these2keep)
        if(length(these2keep) <= 2) next # not enough individuals; skip this batch

        # subset the rest
        K <- subset_kinship(kinship, ind=these2keep)
        ac <- addcovar; if(!is.null(ac)) { ac <- ac[these2keep,,drop=FALSE]; ac <- drop_depcols(ac, TRUE, tol) }
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) { ic <- ic[these2keep,,drop=FALSE]; ic <- drop_depcols(ic, TRUE, tol) }
        wts <- weights; if(!is.null(wts)) wts <- wts[these2keep]
        ph <- pheno[these2keep,phecol,drop=FALSE]

        # multiply stuff by the weights
        K <- weight_kinship(K, wts)
        ac <- weight_matrix(ac, wts)
        Xc <- weight_matrix(Xc, wts)
        ph <- weight_matrix(ph, wts)

        # eigen decomposition of kinship matrix
        Ke <- decomp_kinship(K, cores=cores)

        # fit LMM for each phenotype, one at a time
        if(estimate_hsq) {
            nullresult <- calc_hsq_clean(Ke=Ke, pheno=ph, addcovar=ac, Xcovar=Xc,
                                         is_x_chr=is_x_chr, weights=wts, reml=reml,
                                         cores=cores, check_boundary=check_boundary, tol=tol)

            hsq[, phecol] <- nullresult$hsq
        }
        else {
            # for the log likelihood, calculate the reml=FALSE version
            loglik <- calc_nullLL_clean(hsq=hsq[,phecol,drop=FALSE], Ke=Ke, pheno=ph,
                                        addcovar=ac, Xcovar=Xc,
                                        is_x_chr=is_x_chr, weights=wts, reml=FALSE,
                                        cores=cores)
            nullresult <- list(hsq=hsq[,phecol,drop=FALSE],
                               loglik=loglik)
        }

        # weighted least squares genome scan, using cluster_lapply across chromosomes
        lod <- scan1_pg_clean(genoprobs, these2keep, Ke, ph, ac, ic, is_x_chr,
                              wts, genoprob_Xcol2drop,
                              nullresult$hsq, nullresult$loglik, reml, cores,
                              intcovar_method, tol)

        result[,phecol] <- apply(lod, 2, tapply, chr_index, max, na.rm=TRUE)
    }

    if(!by_chr) result <- apply(result, 2, max, na.rm=TRUE)

    # add attributes
    attr(result, "hsq") <- hsq
    attr(result, "sample_size") <- n

    result
}
