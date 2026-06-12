#' General genome scan with a single-QTL model
#'
#' General genome scan with a single-QTL model, with possible allowance for covariates
#' and a polygeneic effect, but taking a general function that calculates a LOD score
#' at a single position.
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
#' @param func Function to calculate log10 likelihood. It is called as
#' `func(pr, pheno, kinship, addcovar, intcovar, weights, ...)`
#' where `pr` is for the genotype probabilities at a fixed
#' position. It will also be used to calculate the null log10
#' likelihood, with `pr` being NULL in that case.
#' @param vectorize_func If TRUE (the default), assume that `func`
#' takes only a single phenotype. If FALSE, assume `func` is able
#' to take a matrix of phenotypes and return a vector of log10
#' likelihood values.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' @param ... Additional arguments passed to `func()`.
#'
#' @return An object of class `"scan1"`: a matrix of LOD scores, positions x phenotypes.
#'
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
#' # covariates for X chr under null
#' Xcovar <- get_x_covar(iron)
#'
#' # create binary trait
#' bin_pheno <- setNames(as.numeric(iron$pheno[,1] > median(iron$pheno[,1])),
#'                       rownames(iron$pheno))
#'
#' ll_glm <-
#'    function(pr, pheno, addcovar=NULL, ...)
#' {
#'    formula <- ifelse(is.null(pr), "pheno ~ 1", "pheno ~ pr")
#'    if(!is.null(addcovar)) formula <- paste(formula, "+ addcovar")
#'
#'     glm_out <- glm(as.formula(formula), family=binomial(link=probit))
#'     -glm_out$deviance/(2*log(10)) # log10 likelihood
#' }
#'
#' # perform genome scan using glm() with probit link function
#' out <- scan1gen(probs, bin_pheno, Xcovar=Xcovar, func=ll_glm)
#'
#' @seealso [scan1()], [scan1perm()], [scan1max()]
#'
#' @export
scan1gen <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, func, vectorize_func=TRUE,
             cores=1, ...)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    if(!is.null(max_batch) && !is_pos_number(max_batch)) stop("max_batch should be a single positive integer")
    else {
        check_extra_dots(dotargs, c("tol", "quiet", "max_batch"))
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

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
                               weights, complete.cases=TRUE)
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

    phe_batches <- batch_cols(pheno[ind2keep,,drop=FALSE], max_batch)

    # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
    genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # vectorize function if necessary
    if(ncol(pheno)>1 && vectorize_func) {
        vfunc <- function(pr, pheno, ...) {
            if(ncol(pheno) > 1)
                return(sapply(1:ncol(pheno), function(i) func(pr, pheno[,i,drop=FALSE], ...)))
            else return(func(pr, pheno, ...))
        }
    } else {
        vfunc <- func
    }

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

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
        } else pr <- genoprobs[[chr]][these2keep,-1,,drop=FALSE]

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

        # FIX_ME: fitting null model multiple times :(
        ll0 <- vfunc(NULL, ph, addcovar=ac0, weights=wts)

        ll <- matrix(ncol=ncol(ph), nrow=dim(pr)[3])
        for(pos in 1:nrow(ll)) {
            ll[pos,] <- vfunc(pr[,,pos], ph, addcovar=ac, intcovar=ic, weights=wts)
        }

        # calculate LOD score
        lod <- t(ll) - ll0 # deal with multiple phenotypes

        list(lod=lod, n=nrow(ph)) # return LOD & number of individuals used
    }

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- dim(genoprobs)[3,]
    totpos <- sum(npos_by_chr)
    pos_index <- split(seq_len(totpos), rep(seq_len(length(genoprobs)), npos_by_chr))

    # object to contain the LOD scores; also attr to contain sample size
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)
    if(totpos==0) { # edge case of no genoprobs
        colnames(result) <- colnames(pheno)
        attr(result, "sample_size") <- n
        class(result) <- c("scan1", "matrix")
        return(result)
    }

    if(n_cores(cores)==1) { # no parallel processing
        for(i in run_indexes) {
            chr <- run_batches$chr[i]
            chrnam <- names(genoprobs)[chr]
            phebatch <- phe_batches[[run_batches$phe_batch[i]]]
            phecol <- phebatch$cols

            this_result <- by_group_func(i)
            if(!is.null(this_result)) {
                result[pos_index[[chr]], phecol] <- t(this_result$lod)
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
                result[pos_index[[chr]], phecol] <- t(list_result[[i]]$lod)
                if(chr==1) n[phecol] <- list_result[[i]]$n
            }
        }
    }

    pos_names <- unlist(dimnames(genoprobs)[[3]])
    names(pos_names) <- NULL # this is just annoying
    dimnames(result) <- list(pos_names, colnames(pheno))

    # add some attributes with details on analysis
    attr(result, "sample_size") <- n

    class(result) <- c("scan1", "matrix")
    result
}
