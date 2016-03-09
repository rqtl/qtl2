#' Genome scan with a single-QTL and linear mixed model
#'
#' Genome scan with a single-QTL and linear mixed model to account for
#' a random polygenic effect, with possible allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' \code{\link[qtl2geno]{calc_genoprob}}.
#' @param pheno A matrix of phenotypes, individuals x phenotypes.
#' @param kinship A kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional matrix of additive covariates.
#' @param Xcovar An optional matrix with additional additive covariates used for
#' null hypothesis when scanning the X chromosome.
#' @param intcovar An optional matrix of interactive covariates.
#' @param reml If true, use REML; otherwise, use maximimum likelihood.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#' @param ... Additional control parameters; see Details.
#'
#' @return A list containing the following
#' \itemize{
#' \item \code{lod} - A matrix of LOD scores, positions x phenotypes.
#' \item \code{map} - A list containing the map positions at which the
#'     calculations were performed, taken from the input \code{genoprobs}.
#' \item \code{hsq} - A matrix of estimated heritabilities under the
#'     null hypothesis of no QTL. Columns are the phenotypes. If the
#'     \code{"loco"} method was used with
#'     \code{\link[qtl2geno]{calc_kinship}} to calculate a list of kinship
#'     matrices, one per chromosome, the rows of \code{hsq} will be the
#'     heritabilities for the different chromosomes (well, leaving out
#'     each one). If \code{Xcovar} was not NULL, there will at least be an
#'     autosome and X chromosome row.
#' \item \code{addcovar} - Names of additive covariates that were used.
#' \item \code{Xcovar} - Names of special covariates for X chromosome
#'     under the null hypothesis of no QTL
#' \item \code{intcovar} - Names of interactive covariates that were used.
#' \item \code{sample_size} - Vector of sample sizes used for each
#'     phenotype
#' \item \code{snpinfo} - Present only if the input \code{genoprobs}
#'     was produced by \code{\link{genoprob_to_snpprob}}, this is a list
#'     of data frames giving information about all SNPs. The \code{lod}
#'     matrix will contain only results for distinct SNPs. The
#'     \code{index} column in \code{snpinfo} is the row index in the
#'     \code{lod} matrix that corresponds to the current SNP.
#' }
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals. The \code{genoprobs}
#' object should have a component \code{"is_x_chr"} that indicates
#' which of the chromosomes is the X chromosome, if any.
#'
#' If \code{kinship} is a single matrix, then the \code{hsq}
#' attribute in the results is a vector of heritabilities. If
#' \code{kinship} is a list (one matrix per chromosome), then
#' \code{hsq} is a matrix, phenotypes x chromosomes.
#'
#' If \code{reml=TRUE}, restricted maximum likelihood (reml) is used
#' to estimate the heritability under the null hypothesis of no QTL,
#' separately for each phenotype. But then in the genome scans, this
#' is taken as fixed and known, and the usual Gaussian log likelihood
#' is used to calculate LOD scores.
#'
#' The \code{...} argument can contain two additional control
#' parameters; suspended for simplicity (or confusion, depending on
#' your point of view). \code{tol} is used as a tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted);
#' default \code{1e-12}. \code{intcovar_method} indicates whether to
#' use a high-memory (but potentially faster) method or a low-memory
#' (and possibly slower) method, with values \code{"highmem"} or
#' \code{"lowmem"}; default \code{"lowmem"}.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # leave-one-chromosome-out kinship matrices
#' kinship <- calc_kinship(probs, "loco")
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1_lmm(probs, pheno, kinship, covar, Xcovar)
#'
#' @export
scan1_lmm <-
    function(genoprobs, pheno, kinship, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, reml=TRUE, cores=1, ...)
{
    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    stopifnot(tol > 0)
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    check_boundary <- grab_dots(dotargs, "check_boundary", TRUE)
    check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch",
                                "check_boundary"))

    # force things to be matrices
    if(!is.matrix(pheno))
        pheno <- as.matrix(pheno)
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(Xcovar) && !is.matrix(Xcovar))
        Xcovar <- as.matrix(Xcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)

    # check that kinship matrices are square with same IDs
    kinshipIDs <- check_kinship(kinship, length(genoprobs$probs))

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
                               kinshipIDs, complete.cases=TRUE)
    ind2keep <- get_common_ids(ind2keep, rownames(pheno)[rowSums(!is.na(pheno)) > 0])
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
    is_x_chr <- genoprobs$is_x_chr
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs$probs))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- vapply(genoprobs$probs, function(a) dim(a)[3], 1)
    totpos <- sum(npos_by_chr)
    pos_index <- split(1:totpos, rep(seq(along=genoprobs$probs), npos_by_chr))
    pos_names <- unlist(lapply(genoprobs$probs, function(a) dimnames(a)[[3]]))
    names(pos_names) <- NULL # this is just annoying

    # to contain the results
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    dimnames(result) <- list(pos_names, colnames(pheno))

    # number of chr to consider under null
    if(is.list(kinship)) n_null_chr <- length(kinship)
    else if(!is.null(Xcovar)) n_null_chr <- any(is_x_chr) + any(!is_x_chr)
    else n_null_chr <- 1

    hsq <- matrix(nrow=n_null_chr, ncol=ncol(pheno))
    dimnames(hsq) <- hsq_dimnames(kinship, Xcovar, is_x_chr, pheno)
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)

    # loop over batches of phenotypes with the same pattern of NAs
    for(batch in seq(along=phe_batches)) {

        # info about batch
        omit <- phe_batches[[batch]]$omit # ind to omit
        phecol <- phe_batches[[batch]]$cols # phenotype columns in batch

        # individuals to keep in this batch
        these2keep <- ind2keep
        if(length(omit)>0) these2keep <- ind2keep[-omit]
        n[phecol] <- length(these2keep)
        if(length(these2keep) <= 2) next # not enough individuals; skip this batch

        # subset the rest
        if(is.list(kinship)) K <- lapply(kinship, function(a) a[these2keep,these2keep])
        else K <- kinship[these2keep, these2keep]
        ac <- addcovar; if(!is.null(ac)) ac <- ac[these2keep,,drop=FALSE]
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) ic <- ic[these2keep,,drop=FALSE]
        ph <- pheno[these2keep,phecol,drop=FALSE]

        # eigen decomposition of kinship matrix
        Ke <- decomp_kinship(K, cores=cores)

        # fit LMM for each phenotype, one at a time
        nullresult <- calc_hsq_clean(Ke, ph, ac, Xc, is_x_chr, reml, cores,
                                       check_boundary, tol)
        hsq[, phecol] <- nullresult$hsq

        # weighted least squares genome scan, using cluster_lapply across chromosomes
        lod <- scan1_lmm_clean(genoprobs, these2keep, Ke, ph, ac, ic, is_x_chr,
                               genoprob_Xcol2drop,
                               nullresult$hsq, nullresult$loglik, reml, cores,
                               intcovar_method, tol)

        result[,phecol] <- lod
    }

    result <- list(lod = result,
                   map = genoprobs$map,
                   hsq = hsq,
                   sample_size = n,
                   addcovar = colnames4attr(addcovar),
                   Xcovar = colnames4attr(Xcovar),
                   intcovar = colnames4attr(intcovar))

    # preserve any snpinfo from genoprob_to_snpprob
    if("snpinfo" %in% names(genoprobs))
        result$snpinfo <- genoprobs$snpinfo

    class(result) <- c("scan1", "matrix")
    result
}


# fit LMM for each of a matrix of phenotypes
calc_hsq_clean <-
    function(Ke, pheno, addcovar, Xcovar, is_x_chr, reml=TRUE, cores=1,
             check_boundary, tol)
{
    n <- nrow(pheno)
    nphe <- ncol(pheno)

    # if just one kinship matrix, force it to be a list
    if(!is.list(Ke[[1]])) {
        # X chromosome with special covariates
        if(!is.null(Xcovar) && any(is_x_chr) && any(!is_x_chr)) {
            Ke <- list(Ke, Ke)
            is_x_chr <- c(FALSE, TRUE)
        }
        else { Ke <- list(Ke); is_x_chr <- FALSE }
    }

    # function that does the work
    by_chr_func <-
        function(chr)
        {
            # premultiply phenotypes and covariates by transposed eigenvectors
            y <- Ke[[chr]]$vectors %*% pheno
            ac <- cbind(rep(1, n), addcovar)
            if(!is.null(Xcovar) && is_x_chr[chr]) # add Xcovar if necessary
                ac <- drop_depcols(cbind(ac, Xcovar), FALSE, tol)
            logdetXpX = Rcpp_calc_logdetXpX(ac)
            ac <- Ke[[chr]]$vectors %*% ac

            Rcpp_fitLMM_mat(Ke[[chr]]$values, y, ac, reml, check_boundary,
                            logdetXpX, tol)
        }

    # now do the work
    result <- cluster_lapply(cores, seq(along=Ke), by_chr_func)

    # re-arrange results
    hsq <- matrix(unlist(lapply(result, function(a) a$hsq)), byrow=TRUE,
                  ncol=nphe)
    loglik <- matrix(unlist(lapply(result, function(a) a$loglik)), byrow=TRUE,
                  ncol=nphe)

    list(hsq=hsq, loglik=loglik)
}

# perform the LMM scan
# genoprobs is still a big complicated calc_genoprob object
scan1_lmm_clean <-
    function(genoprobs, ind2keep, Ke, pheno, addcovar, intcovar, is_x_chr,
             genoprob_Xcol2drop,
             hsq, null_loglik, reml, cores, intcovar_method, tol)
{
    n <- nrow(pheno)
    nphe <- ncol(pheno)

    if(!is.list(Ke[[1]])) {
        loco <- FALSE
        if(nrow(hsq)==1) no_x <- TRUE
        else no_x <- FALSE
    } else loco <- TRUE

    batches <- list(chr=rep(seq(along=genoprobs$probs), ncol(pheno)),
                    phecol=rep(1:ncol(pheno), each=length(genoprobs$probs)))

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
            ac <- cbind(rep(1, n), addcovar)
            ic <- intcovar

            # subset the genotype probabilities: drop cols with all 0s, plus the first column
            Xcol2drop <- genoprob_Xcol2drop[[chr]]
            if(length(Xcol2drop) > 0) {
                pr <- genoprobs$probs[[chr]][ind2keep,-Xcol2drop,,drop=FALSE]
                pr <- pr[,-1,,drop=FALSE]
            }
            else
                pr <- genoprobs$probs[[chr]][ind2keep,-1,,drop=FALSE]

            # calculate weights for this chromosome
            if(loco) {
                weights <- 1/(hsq[chr,phecol]*Keval + (1-hsq[chr,phecol]))
                nullLL <- null_loglik[chr,phecol]
            }
            else {
                if(no_x || !is_x_chr[chr]) {
                    weights <- 1/(hsq[1,phecol]*Keval + (1-hsq[1,phecol]))
                    nullLL <- null_loglik[1,phecol]
                }
                else {
                    weights <- 1/(hsq[2,phecol]*Keval + (1-hsq[2,phecol]))
                    nullLL <- null_loglik[2,phecol]
                }
            }
            weights <- sqrt(weights)

            if(is.null(ic))
                loglik <- scan_lmm_onechr(pr, y, ac, Kevec, weights, tol)
            else if(intcovar_method=="highmem")
                loglik <- scan_lmm_onechr_intcovar_highmem(pr, y, ac, ic, Kevec, weights, tol)
            else
                loglik <- scan_lmm_onechr_intcovar_lowmem(pr, y, ac, ic, Kevec, weights, tol)
            lod <- (loglik - nullLL)/log(10)
        }

    # now do the work
    lod_list <- cluster_lapply(cores, seq(along=batches$chr), by_batch_func)

    npos_by_chr <- vapply(genoprobs$probs, function(a) dim(a)[3], 1)
    totpos <- sum(npos_by_chr)
    pos_index <- split(1:totpos, rep(seq(along=genoprobs$probs), npos_by_chr))

    # to contain the results
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    for(batch in seq(along=batches$chr)) {
        chr <- batches$chr[batch]
        phecol <- batches$phecol[batch]
        result[pos_index[[chr]], phecol] <- lod_list[[batch]]
    }

    result
}

# check that kinship matrices are square
# and have same row and column IDs
#
# returns vector of IDs
check_kinship <-
    function(kinship, n_chr)
{
    if(!is.list(kinship)) { # one kinship matrix
        stopifnot(nrow(kinship) == ncol(kinship))
        stopifnot( all(rownames(kinship) == colnames(kinship)) )
        return(rownames(kinship))
    } else {
        if(length(kinship) != n_chr)
            stop("length(kinship) != no. chromosomes (", n_chr, ")")

        kinship_square <- vapply(kinship, function(mat) nrow(mat) == ncol(mat), TRUE)
        stopifnot( all(kinship_square) )

        kinship_sameIDs <- vapply(kinship, function(mat) (nrow(mat) == nrow(kinship[[1]])) &&
                                  all((rownames(mat) == rownames(kinship[[1]])) &
                                      (colnames(mat) == colnames(kinship[[1]])) &
                                      (rownames(mat) == colnames(kinship[[1]]))), TRUE)
        if(!all(kinship_sameIDs))
            stop("All kinship matrices should be the same size ",
                 "and have the same row and column names")

        return(rownames(kinship))
    }
}

# dimnames for hsq
# (a bit awkward due to loco or not, and xchr special or not)
hsq_dimnames <-
    function(kinship, Xcovar, is_x_chr, pheno)
{
    if(is.list(kinship)) rn <- names(kinship)
    else if(!is.null(Xcovar) && any(is_x_chr) + any(!is_x_chr) == 2) {
        rn <- c(names(is_x_chr)[!is_x_chr][1],
                names(is_x_chr)[is_x_chr][1])
    }
    else {
        rn <- names(is_x_chr)[!is_x_chr][1]
    }

    list(rn, colnames(pheno))
}
