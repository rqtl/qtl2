#' Genome scan with a single-QTL and linear mixed model
#'
#' Genome scan with a single-QTL and linear mixed model to account for
#' a random polygenic effect, with possible allowance for covariates.
#'
#' @param genoprobs A list of 3-dimensional arrays of genotype
#' probabilities; each component is a chromosome, and has dimension
#' individuals x genotypes x positions.
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
#' @return A matrix of LOD scores, positions x phenotypes.
#' Heritabilities (estimated under the null hypothesis, of no QTL) are
#' included as an attribute \code{herit}. Covariate column names are
#' included as attributes (\code{"addcovar"}, \code{"intcovar"}, and
#' \code{"Xcovar"}), as is a vector with the sample size for each
#' phenotype (\code{"sample_size"})
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals. The \code{genoprobs}
#' object should have an attribute \code{"is_x_chr"} that indicates
#' which of the chromosomes is the X chromosome, if any.
#'
#' If \code{kinship} is a single matrix, then the \code{herit}
#' attribute in the results is a vector of heritabilities. If
#' \code{kinship} is a list (one matrix per chromosome), then
#' \code{herit} is a matrix, phenotypes x chromosomes.
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
    if("tol" %in% names(dotargs))
        tol <- dotargs$tol
    else tol <- 1e-12

    if("intcovar_method" %in% names(dotargs)) {
        intcovar_method <- dotargs$intcovar_method
        if(!(intcovar_method %in% c("highmem", "lowmem"))) {
            warning('intcovar_method "', intcovar_method, '" not valid; using "lowmem".')
            intcovar_method <- "lowmem"
        }
    } else intcovar_method <- "lowmem"

    if("quiet" %in% names(dotargs))
        quiet <- dotargs$quiet
    else quiet <- TRUE

    if("max_batch" %in% names(dotargs))
        max_batch <- dotargs$max_batch
    else max_batch <- NULL

    if("check_boundary" %in% names(dotargs))
        check_boundary <- dotargs$check_boundary
    else check_boundary <- TRUE

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
    kinshipIDs <- check_kinship(kinship)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs[[1]], addcovar, Xcovar, intcovar,
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
    class(genoprobs) <- "list" # treat as regular list
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- vapply(genoprobs, function(a) dim(a)[3], 1)
    totpos <- sum(npos_by_chr)
    pos_index <- split(1:totpos, rep(seq(along=genoprobs), npos_by_chr))
    pos_names <- unlist(lapply(genoprobs, function(a) dimnames(a)[[3]]))

    # to contain the results
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    dimnames(result) <- list(pos_names, colnames(pheno))

    # number of chr to consider under null
    if(is.list(kinship)) {
        n_null_chr <- length(kinship)
        null_chrnames <- names(kinship)
    }
    else if(!is.null(Xcovar) && any(is_x_chr) && any(!is_x_chr)) {
        n_null_chr <- 2
        null_chrnames <- c(names(is_x_chr)[!is_x_chr][1],
                           names(is_x_chr)[is_x_chr][1])
    }
    else {
        n_null_chr <- 1
        null_chrnames <- names(genoprobs)[1]
    }

    hsq <- matrix(nrow=n_null_chr, ncol=ncol(pheno))
    dimnames(hsq) <- list(null_chrnames, colnames(pheno))
    null_loglik <- hsq
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
        nullresult <- calc_hsq_clean(Ke, ph, ac, Xc, is_x_chr, reml, cores=cores,
                                       check_boundary, tol)
        hsq[, phecol] <- nullresult$hsq
        null_loglik[, phecol] <- nullresult$loglik

        # weighted least squares genome scan, using cluster_lapply across chromosomes

    }

    # temporary result
    result <- list(hsq=hsq, null_loglik=null_loglik)

    attr(result, "sample_size") <- n
    attr(result, "addcovar") <- colnames4attr(addcovar)
    attr(result, "Xcovar") <- colnames4attr(Xcovar)
    attr(result, "intcovar") <- colnames4attr(intcovar)

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
        if(any(is_x_chr) && !is.null(Xcovar)) { # X chromosome with special covariates
            Ke <- list(Ke, Ke)
            names(Ke) <- c(names(is_x_chr)[!is_x_chr][1],
                           names(is_x_chr)[is_x_chr][1])
        }
        else {
            Ke <- list(Ke)
            names(Ke) <- names(is_x_chr)[1]
        }
    }

    # function that does the work
    by_chr_func <-
        function(chr)
        {
            # premultiply phenotypes and covariates by transposed eigenvectors
            y <- Ke[[chr]]$vectors %*% pheno
            ac <- cbind(rep(1, n), addcovar)
            if(!is.null(Xcovar) && is_x_chr[names(Ke)[chr]]) { # add Xcovar if necessary
                ac <- drop_depcols(cbind(ac, Xcovar), FALSE, tol)
            }
            ac <- Ke[[chr]]$vectors %*% ac
            logdetXpX = Rcpp_calc_logdetXpX(ac)

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


# check that kinship matrices are square
# and have same row and column IDs
#
# returns vector of IDs
check_kinship <-
    function(kinship)
{
    if(!is.list(kinship)) { # one kinship matrix
        stopifnot(nrow(kinship) == ncol(kinship))
        stopifnot( all(rownames(kinship) == colnames(kinship)) )
        return(rownames(kinship))
    } else {
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
