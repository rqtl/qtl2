#' Maximum LOD score from genome scan with a single-QTL model
#'
#' Maximum LOD score from genome scan with a single-QTL model by
#' Haley-Knott regression or a linear mixed model, with possible
#' allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param pheno A numeric matrix of phenotypes, individuals x phenotypes.
#' @param kinship Optional kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional numeric matrix of additive covariates.
#' @param Xcovar An optional numeric matrix with additional additive covariates used for
#' null hypothesis when scanning the X chromosome.
#' @param intcovar An numeric optional matrix of interactive covariates.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param reml If `kinship` provided: if `reml=TRUE`, use
#' REML; otherwise maximum likelihood.
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param hsq Considered only if `kinship` is provided, in which case
#'     this is taken as the assumed value for the residual
#'     heritability. It should be a vector with length corresponding
#'     to the number of columns in `pheno`, or (if `kinship`
#'     corresponds to a list of LOCO kinship matrices) a matrix with dimension
#'     `length(kinship) x ncol(pheno)`.
#' @param by_chr If TRUE, save the individual chromosome maxima.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#' @param ... Additional control parameters; see Details.
#'
#' @return Either a vector of genome-wide maximum LOD scores, or if
#'     `by_chr` is TRUE, a matrix with the chromosome-specific maxima,
#'     with the rows being the chromosomes and the columns being the
#'     phenotypes.
#'
#' @details Equivalent to running `scan1()` and then saving the column
#'     maxima, with some savings in memory usage.
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c("19", "X")] # subset to chr 19 and X}
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1max(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' @seealso [scan1()], [scan1perm()]
#'
#' @export
scan1max <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE,
             model=c("normal", "binary"), hsq=NULL,
             by_chr=FALSE, cores=1, ...)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    # grab dot args
    dotargs <- list(...)
    if("n_perm" %in% names(dotargs))
        stop("You included n_perm as an argument; you probably want to run scan1perm not scan1.")

    model <- match.arg(model)

    if(!is.null(kinship)) { # fit linear mixed model
        if(model=="binary") warning("Can't fit binary model with kinship matrix; using normal model")
        return(scan1max_pg(genoprobs, pheno, kinship, addcovar, Xcovar, intcovar,
                           weights, reml, hsq, by_chr, cores, ...))
    }

    # deal with the dot args
    tol <- grab_dots(dotargs, "tol", 1e-12)
    if(!is_pos_number(tol)) stop("tol should be a single positive number")
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    if(!is.null(max_batch) && !is_pos_number(max_batch)) stop("max_batch should be a single positive integer")
    if(model=="binary") {
        bintol <- grab_dots(dotargs, "bintol", sqrt(tol)) # for model="binary"
        if(!is_pos_number(bintol)) stop("bintol should be a single positive number")
        eta_max <- grab_dots(dotargs, "eta_max", log(1-tol)-log(tol)) # for model="binary"
        if(!is_pos_number(eta_max)) stop("eta_max should be a single positive number")
        maxit <- grab_dots(dotargs, "maxit", 100) # for model="binary"
        if(!is_nonneg_number(maxit)) stop("maxit should be a single non-negative integer")
        check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch", "maxit", "bintol", "eta_max"))
    }
    else {
        check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch"))
    }

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

    # for binary model
    if(model=="binary") {
        if(!is.null(kinship))
            stop("Can't yet account for kinship with model = \"binary\"")
        pheno <- check_binary_pheno(pheno)
    }
    else {
        # square-root of weights (only if model="normal")
        weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
                               weights, complete.cases=TRUE)
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

    n_chr <- length(genoprobs)

    # batches for analysis, to allow parallel analysis
    run_batches <- data.frame(chr=rep(seq_len(length(genoprobs)), length(phe_batches)),
                              phe_batch=rep(seq_along(phe_batches), each=length(genoprobs)))
    run_indexes <- seq_len(length(genoprobs)*length(phe_batches))

    # the function that does the work
    by_group_func <- function(i) {
        # deal with batch information, including individuals to drop due to missing phenotypes
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        phecol <- phebatch$cols
        omit <- phebatch$omit
        these2keep <- ind2keep # individuals 2 keep for this batch
        if(length(omit) > 0) these2keep <- ind2keep[-omit]
        if(length(these2keep)<=2) return(NULL) # not enough individuals

        # subset the genotype probabilities: drop cols with all 0s, plus the first column
        Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
        if(length(Xcol2drop) > 0) {
            pr <- genoprobs[[chr]][these2keep,-Xcol2drop,,drop=FALSE]
            pr <- pr[,-1,,drop=FALSE]
        }
        else
            pr <- genoprobs[[chr]][these2keep,-1,,drop=FALSE]

        # subset the rest
        ac <- addcovar; if(!is.null(ac)) { ac <- ac[these2keep,,drop=FALSE]; ac <- drop_depcols(ac, TRUE, tol) }
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) { ic <- ic[these2keep,,drop=FALSE]; ic <- drop_depcols(ic, TRUE, tol) }
        ph <- pheno[these2keep,phecol,drop=FALSE]
        wts <- weights[these2keep]

        # if X chr, paste X covariates onto additive covariates
        # (only for the null)
        if(is_x_chr[chr]) ac0 <- drop_depcols(cbind(ac, Xc), add_intercept=FALSE, tol)
        else ac0 <- ac

        if(model=="normal") {
            # FIX_ME: calculating null RSS multiple times :(
            nullrss <- nullrss_clean(ph, ac0, wts, add_intercept=TRUE, tol)

            # scan1 function taking clean data (with no missing values)
            rss <- scan1_clean(pr, ph, ac, ic, wts, add_intercept=TRUE, tol, intcovar_method)

            # calculate LOD score
            lod <- nrow(ph)/2 * (log10(nullrss) - log10(rss))
        }
        else { # binary traits
            # FIX_ME: calculating null LOD multiple times :(
            nulllod <- null_binary_clean(ph, ac0, wts, add_intercept=TRUE, maxit, bintol, tol, eta_max)

            # scan1 function taking clean data (with no missing values)
            lod <- scan1_binary_clean(pr, ph, ac, ic, wts, add_intercept=TRUE,
                                      maxit, bintol, tol, intcovar_method, eta_max)

            # calculate LOD score
            lod <- lod - nulllod
        }

        list(lod=apply(lod, 1, max, na.rm=TRUE), n=nrow(ph)) # return LOD & number of individuals used
    }

    # object to contain the LOD scores; also attr to contain sample size
    result <- matrix(nrow=n_chr, ncol=ncol(pheno))
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)
    if(n_chr==0) { # edge case of no genoprobs
        return(NULL)
    }

    if(n_cores(cores)==1) { # no parallel processing
        for(i in run_indexes) {
            chr <- run_batches$chr[i]
            chrnam <- names(genoprobs)[chr]
            phebatch <- phe_batches[[run_batches$phe_batch[i]]]
            phecol <- phebatch$cols

            this_result <- by_group_func(i)
            if(!is.null(this_result)) {
                result[chr, phecol] <- this_result$lod
                if(chr==1) n[phecol] <- this_result$n
            }
        }
    }
    else {
        # calculations in parallel
        list_result <- cluster_lapply(cores, run_indexes, by_group_func)

        # check for problems (if clusters run out of memory, they'll return NULL)
        result_is_null <- vapply(list_result, is.null, TRUE)
        if(any(result_is_null))
            stop("cluster problem: returned ", sum(result_is_null), " NULLs.")

        # reorganize results
        for(i in run_indexes) {
            chr <- run_batches$chr[i]
            chrnam <- names(genoprobs)[chr]
            phebatch <- phe_batches[[run_batches$phe_batch[i]]]
            phecol <- phebatch$cols

            if(!is.null(list_result[[i]])) {
                result[chr, phecol] <- list_result[[i]]$lod
                if(chr==1) n[phecol] <- list_result[[i]]$n
            }
        }
    }

    if(by_chr) {
        dimnames(result) <- list(names(genoprobs), colnames(pheno))
    } else {
        result <- apply(result, 2, max, na.rm=TRUE)
        names(result) <- colnames(pheno)
    }

    # add some attributes with details on analysis
    attr(result, "sample_size") <- n

    result
}
