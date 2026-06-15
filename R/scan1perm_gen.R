# scan1perm with general function but no covariates, kinship, weights
scan1perm_gen_simple <-
    function(genoprobs, pheno, n_perm=1, perm_strata=NULL,
             cores=1, scan_func, ind2keep, ..., max_batch=NULL)
{
    dotargs <- list(...)
    if("quiet" %in% names(dotargs)) quiet <- dotargs$quiet
    else quiet <- TRUE

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # max batch size
    if(is.null(max_batch))
        max_batch <- min(1000, ceiling(n_perm*length(genoprobs)*ncol(pheno)/n_cores(cores)))

    # generate permutations
    perms <- gen_strat_perm(n_perm, ind2keep, perm_strata)

    # batch permutations
    phe_batches <- batch_vec( rep(seq_len(ncol(pheno)), n_perm), max_batch)
    perm_batches <- batch_vec( rep(seq_len(n_perm), each=ncol(pheno)), max_batch)

    # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
    genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # batches for analysis, to allow parallel analysis
    run_batches <- data.frame(chr=rep(seq_len(length(genoprobs)), length(phe_batches)),
                              phe_batch=rep(seq_along(phe_batches), each=length(genoprobs)))

    run_indexes <- seq_len(length(genoprobs)*length(phe_batches))

    # subset the phenotypes
    pheno <- pheno[ind2keep,,drop=FALSE]
    nind <- nrow(pheno)

    # the function that does the work
    by_group_func <- function(i) {
        # deal with batch information, including individuals to drop due to missing phenotypes
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        permbatch <- perm_batches[[run_batches$phe_batch[i]]]

        # subset the genotype probabilities: drop cols with all 0s, plus the first column
        Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
        if(length(Xcol2drop) > 0)
            pr <- genoprobs[[chr]][ind2keep,-Xcol2drop,,drop=FALSE]
        else
            pr <- genoprobs[[chr]][ind2keep,-1,,drop=FALSE]
        pr <- list("chr"=pr)
        class(this_pr) <- class(genoprobs)

        ph <- pheno[,phebatch,drop=FALSE]
        for(col in seq_len(ncol(ph))) # permute columns
            ph[,col] <- ph[perms[,permbatch[col]] , col]

        lod <- scan_func(genoprobs=pr, pheno=ph,
                         addcovar=NULL, Xcovar=NULL,
                         intcovar=NULL, kinship=NULL,
                         cores=1, ...)

        # return column maxima
        apply(lod, 2, max, na.rm=TRUE)
    }

    # calculations in parallel
    list_result <- cluster_lapply(cores, run_indexes, by_group_func)

    # check for problems (if clusters run out of memory, they'll return NULL)
    result_is_null <- vapply(list_result, is.null, TRUE)
    if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

    # reorganize results
    result <- array(dim=c(length(genoprobs), n_perm, ncol(pheno)))

    for(i in run_indexes) {
        chr <- run_batches$chr[i]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        permbatch <- perm_batches[[run_batches$phe_batch[i]]]

        for(j in seq_along(phebatch))
            result[chr,permbatch[j], phebatch[j]] <- list_result[[i]][j]

    }

    result <- apply(result, c(2,3), max)
    colnames(result) <- colnames(pheno)

    class(result) <- c("scan1perm", "matrix")
    result

}

scan1perm_gen <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL,
             n_perm=1, perm_strata=NULL,
             cores=1, scan_func, ind2keep, ..., max_batch=NULL)
{
    dotargs <- list(...)
    if("quiet" %in% names(dotargs)) quiet <- dotargs$quiet
    else quiet <- TRUE

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # max batch size
    if(is.null(max_batch))
        max_batch <- min(1000, ceiling(n_perm*length(genoprobs)*ncol(pheno)/n_cores(cores)))

    # generate permutations
    perms <- gen_strat_perm(n_perm, ind2keep, perm_strata)

    # batch permutations
    phe_batches <- batch_cols(pheno[ind2keep,,drop=FALSE], max_batch)

    # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
    genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # batches for analysis, to allow parallel analysis
    run_batches <- data.frame(chr=rep(seq_len(length(genoprobs)), length(phe_batches)*n_perm),
                              phe_batch=rep(seq_along(phe_batches), each=length(genoprobs)*n_perm),
                              perm_batch=rep(rep(seq_len(n_perm), each=length(genoprobs), length(phe_batches))))

    run_indexes <- seq_len(length(genoprobs)*length(phe_batches)*n_perm)

    # the function that does the work
    by_group_func <- function(i) {
        # deal with batch information, including individuals to drop due to missing phenotypes
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        permbatch <- run_batches$perm_batch[i]
        phecol <- phebatch$cols
        omit <- phebatch$omit
        these2keep <- ind2keep # individuals 2 keep for this batch
        if(length(omit) > 0) these2keep <- ind2keep[-omit]
        if(length(these2keep)<=2) return(NULL) # not enough individuals

        # apply permutation to the probs
        #     (I initially thought we'd need to apply these just to the rownames, but
        #     actually we've already aligned probs and pheno and so forth
        #     via ind2keep/these2keep, so now we need to just permute the rows)
        #
        # (also need some contortions here to keep the sizes the same)
        pr <- genoprobs[[chr]][ind2keep,,,drop=FALSE]
        pr <- pr[perms[,permbatch],,,drop=FALSE]
        rownames(pr) <- ind2keep
        pr <- pr[these2keep,,,drop=FALSE]

        # subset the genotype probabilities: drop cols with all 0s, plus the first column
        Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
        if(length(Xcol2drop) > 0) {
            pr <- pr[,-Xcol2drop,,drop=FALSE]
        }
        pr <- list("chr"=pr)
        class(this_pr) <- class(genoprobs)

        # subset the rest
        ac <- addcovar; if(!is.null(ac)) ac <- ac[these2keep,,drop=FALSE]
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) ic <- ic[these2keep,,drop=FALSE]
        wts <- weights; if(!is.null(wts)) wts <- wts[these2keep]
        k <- kinship; if(is_kinship_list(k)) k <- k[[chr]]
        if(!is.null(k)) k <- k[these2keep,these2keep]
        ph <- pheno[these2keep, phecol, drop=FALSE]

        # if X chr, paste X covariates onto additive covariates
        # (only for the null)
        if(is_x_chr[chr]) ac0 <- drop_depcols(cbind(ac, Xc), add_intercept=FALSE, tol)
        else ac0 <- ac

        # scan1 function taking clean data (with no missing values)
        lod <- scan_func(genoprobs=pr, pheno=ph, kinship=k, addcovar=ac,
                         intcovar=ic, weights=wts, ...)

        # return column maxima
        apply(lod, 2, max, na.rm=TRUE)
    }

    # calculations in parallel
    list_result <- cluster_lapply(cores, run_indexes, by_group_func)

    # check for problems (if clusters run out of memory, they'll return NULL)
    result_is_null <- vapply(list_result, is.null, TRUE)
    if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

    # reorganize results
    result <- array(dim=c(length(genoprobs), n_perm, ncol(pheno)))

    for(i in run_indexes) {
        chr <- run_batches$chr[i]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        permbatch <- run_batches$perm_batch[i]

        result[chr, permbatch, phebatch$cols] <- list_result[[i]]
    }

    result <- apply(result, c(2,3), max)
    colnames(result) <- colnames(pheno)

    class(result) <- c("scan1perm", "matrix")
    result
}



# generate (potentially) stratified permutations
#   (actual work is done in c++, see src/random.cpp)
gen_strat_perm <-
    function(n_perm, ind2keep, perm_strata=NULL)
{
    if(!is.null(perm_strata)) {
        perm_strata <- perm_strata[ind2keep]
        u <- unique(perm_strata)

        if(length(u) > 1) {

            strat_numeric <- match(perm_strata, u)-1

            return(permute_nvector_stratified(n_perm,
                                              seq_along(ind2keep),
                                              strat_numeric,
                                              length(u)))
        }
    }

    permute_nvector(n_perm, seq_along(ind2keep))
}
