# Genome scan with a single-QTL and linear mixed model
#
# called by scan1() when kinship() is provided.
scan1_pg <-
    function(genoprobs, pheno, kinship, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, hsq=NULL,
             cores=1, ...)
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
    ind2keep <- get_common_ids(ind2keep, rownames(pheno)[rowSums(is.finite(pheno)) > 0])
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

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- dim(genoprobs)[3,]
    totpos <- sum(npos_by_chr)
    pos_index <- split(seq_len(totpos), rep(seq_len(length(genoprobs)), npos_by_chr))
    pos_names <- unlist(dimnames(genoprobs)[[3]])
    names(pos_names) <- NULL # this is just annoying

    # to contain the results
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    dimnames(result) <- list(pos_names, colnames(pheno))

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
            loglik <- calc_nullLL_clean(Ke=Ke, pheno=ph, addcovar=ac, Xcovar=Xc,
                                        is_x_chr=is_x_chr, weights=wts, reml=reml,
                                        hsq=hsq[,phecol], cores=cores)
            nullresult <- list(hsq=hsq[,phecol],
                               loglik=loglik)
        }

        # weighted least squares genome scan, using cluster_lapply across chromosomes
        lod <- scan1_pg_clean(genoprobs, these2keep, Ke, ph, ac, ic, is_x_chr,
                              wts, genoprob_Xcol2drop,
                              nullresult$hsq, nullresult$loglik, reml, cores,
                              intcovar_method, tol)

        result[,phecol] <- lod
    }

    # add attributes
    attr(result, "hsq") <- hsq
    attr(result, "sample_size") <- n

    class(result) <- c("scan1", "matrix")
    result
}


# fit LMM for each of a matrix of phenotypes
# Ke is eigendecomposition of 2*kinship
calc_hsq_clean <-
    function(Ke, pheno, addcovar=NULL, Xcovar=NULL, is_x_chr=FALSE, weights=NULL,
             reml=TRUE, cores=1, check_boundary=FALSE, tol=1e-12)
{
    n <- nrow(pheno)
    nphe <- ncol(pheno)

    # if just one kinship matrix, force it to be a list
    if(!is_kinship_list(Ke)) {
        # X chromosome with special covariates
        if(!is.null(Xcovar) && any(is_x_chr) && any(!is_x_chr)) {
            Ke <- list(Ke, Ke)
            is_x_chr <- c(FALSE, TRUE)
        }
        else { Ke <- list(Ke) }
    }

    # function that does the work
    by_chr_func <-
        function(chr)
        {
            # premultiply phenotypes and covariates by transposed eigenvectors
            y <- Ke[[chr]]$vectors %*% pheno
            intercept <- weights; if(is_null_weights(weights)) intercept <- rep(1,n)
            ac <- cbind(intercept, addcovar)
            if(!is.null(Xcovar) && is_x_chr[chr]) # add Xcovar if necessary
                ac <- drop_depcols(cbind(ac, Xcovar), FALSE, tol)
            logdetXpX = Rcpp_calc_logdetXpX(ac)
            ac <- Ke[[chr]]$vectors %*% ac

            Rcpp_fitLMM_mat(Ke[[chr]]$values, y, ac, reml, check_boundary,
                            logdetXpX, tol)
        }

    # now do the work
    result <- cluster_lapply(cores, seq_along(Ke), by_chr_func)

    # check for problems (if clusters run out of memory, they'll return NULL)
    result_is_null <- vapply(result, is.null, TRUE)
    if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

    # re-arrange results
    hsq <- matrix(unlist(lapply(result, function(a) a$hsq)), byrow=TRUE,
                  ncol=nphe)
    loglik <- matrix(unlist(lapply(result, function(a) a$loglik)), byrow=TRUE,
                  ncol=nphe)

    list(hsq=hsq, loglik=loglik)
}

# calculate log likelihood for given hsq and
# for each column of a matrix of phenotypes
# Ke is eigendecomposition of 2*kinship
#
# hsq is a matrix with ncol = ncol(pheno)
calc_nullLL_clean <-
    function(Ke, pheno, addcovar=NULL, Xcovar=NULL, is_x_chr=FALSE, weights=NULL,
             reml=TRUE, hsq, cores=1)
{
    n <- nrow(pheno)
    nphe <- ncol(pheno)
    if(ncol(hsq) != nphe) stop("ncol(hsq) != ncol(pheno)")

    # if just one kinship matrix, force it to be a list
    if(!is_kinship_list(Ke)) {
        # X chromosome with special covariates
        if(!is.null(Xcovar) && any(is_x_chr) && any(!is_x_chr)) {
            Ke <- list(Ke, Ke)
            is_x_chr <- c(FALSE, TRUE)
        }
        else { Ke <- list(Ke) }
    }
    if(nrow(hsq) != length(Ke)) stop("nrow(hsq) != no. chr")

    # function that does the work
    by_chr_func <-
        function(chr)
        {
            # premultiply phenotypes and covariates by transposed eigenvectors
            y <- Ke[[chr]]$vectors %*% pheno
            intercept <- weights; if(is_null_weights(weights)) intercept <- rep(1,n)
            ac <- cbind(intercept, addcovar)
            if(!is.null(Xcovar) && is_x_chr[chr]) # add Xcovar if necessary
                ac <- drop_depcols(cbind(ac, Xcovar), FALSE, tol)
            logdetXpX = Rcpp_calc_logdetXpX(ac)
            ac <- Ke[[chr]]$vectors %*% ac

            Rcpp_calcLL_mat(Ke[[chr]]$values, y, ac, reml, hsq[chr,], logdetXpX)
        }

    # now do the work
    result <- cluster_lapply(cores, seq_along(Ke), by_chr_func)

    # check for problems (if clusters run out of memory, they'll return NULL)
    result_is_null <- vapply(result, is.null, TRUE)
    if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

    # re-arrange results
    matrix(unlist(lapply(result, function(a) a$loglik)), byrow=TRUE,
           ncol=nphe)
}

# perform the LMM scan
# genoprobs is still a big complicated calc_genoprob object
scan1_pg_clean <-
    function(genoprobs, ind2keep, Ke, pheno, addcovar, intcovar, is_x_chr,
             weights, genoprob_Xcol2drop,
             hsq, null_loglik, reml, cores, intcovar_method, tol)
{
    n <- nrow(pheno)
    nphe <- ncol(pheno)

    if(!is.list(Ke[[1]])) {
        loco <- FALSE
        if(nrow(hsq)==1) no_x <- TRUE
        else no_x <- FALSE
    } else loco <- TRUE

    batches <- list(chr=rep(seq_len(length(genoprobs)), ncol(pheno)),
                    phecol=rep(seq_len(ncol(pheno)), each=length(genoprobs)))

    # function that does the work
    by_batch_func <-
        function(batch)
        {
            chr <- batches$chr[batch]
            phecol <- batches$phecol[batch]

            if(loco) {
                Kevec <- Ke[[chr]]$vectors
                Keval <- Ke[[chr]]$values
            }
            else {
                Kevec <- Ke$vectors
                Keval <- Ke$values
            }

            # prep phenotype and covariates
            y <- pheno[,phecol,drop=FALSE]
            intercept <- weights; if(is_null_weights(weights)) intercept <- rep(1,n)
            ac <- cbind(intercept, addcovar)
            ic <- intcovar

            # subset the genotype probabilities: drop cols with all 0s, plus the first column
            Xcol2drop <- genoprob_Xcol2drop[[chr]]
            if(length(Xcol2drop) > 0) {
                pr <- genoprobs[[chr]][ind2keep,-Xcol2drop,,drop=FALSE]
                pr <- pr[,-1,,drop=FALSE]
            }
            else
                pr <- genoprobs[[chr]][ind2keep,-1,,drop=FALSE]
            # weight the probabilities
            pr <- weight_array(pr, weights)

            # calculate weights for this chromosome
            if(loco) {
                lmm_wts <- 1/(hsq[chr,phecol]*Keval + (1-hsq[chr,phecol]))
                nullLL <- null_loglik[chr,phecol]
            }
            else {
                if(no_x || !is_x_chr[chr]) {
                    lmm_wts <- 1/(hsq[1,phecol]*Keval + (1-hsq[1,phecol]))
                    nullLL <- null_loglik[1,phecol]
                }
                else {
                    lmm_wts <- 1/(hsq[2,phecol]*Keval + (1-hsq[2,phecol]))
                    nullLL <- null_loglik[2,phecol]
                }
            }
            lmm_wts <- sqrt(lmm_wts)

            if(is.null(ic))
                loglik <- scan_pg_onechr(pr, y, ac, Kevec, lmm_wts, tol)
            else if(intcovar_method=="highmem")
                loglik <- scan_pg_onechr_intcovar_highmem(pr, y, ac, ic, Kevec, lmm_wts, tol)
            else
                loglik <- scan_pg_onechr_intcovar_lowmem(pr, y, ac, ic, Kevec, lmm_wts, tol)
            lod <- (loglik - nullLL)/log(10)
        }

    # now do the work
    lod_list <- cluster_lapply(cores, seq_along(batches$chr), by_batch_func)

    # check for problems (if clusters run out of memory, they'll return NULL)
    result_is_null <- vapply(lod_list, is.null, TRUE)
    if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

    npos_by_chr <- dim(genoprobs)[3,]
    totpos <- sum(npos_by_chr)
    pos_index <- split(seq_len(totpos), rep(seq_len(length(genoprobs)), npos_by_chr))

    # to contain the results
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    for(batch in seq_along(batches$chr)) {
        chr <- batches$chr[batch]
        phecol <- batches$phecol[batch]
        result[pos_index[[chr]], phecol] <- lod_list[[batch]]
    }

    result
}


# dimnames for hsq
# (a bit awkward due to loco or not, and xchr special or not)
hsq_dimnames <-
    function(kinship, Xcovar, is_x_chr, pheno)
{
    if(is_kinship_list(kinship)) {
        rn <- names(kinship)
    }
    else if(!is.null(Xcovar) && any(is_x_chr) + any(!is_x_chr) == 2) {
        rn <- c(names(is_x_chr)[!is_x_chr][1],
                names(is_x_chr)[is_x_chr][1])
    }
    else {
        rn <- names(is_x_chr)[!is_x_chr][1]
    }

    list(rn, colnames(pheno))
}
