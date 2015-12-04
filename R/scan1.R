#' Genome scan with a single-QTL model by Haley-Knott regression
#'
#' Genome scan with a single-QTL model by Haley-Knott regression, with possible allowance for covariates.
#'
#' @param genoprobs A list of 3-dimensional arrays of genotype
#' probabilities; each component is a chromosome, and has dimension
#' individuals x genotypes x positions.
#' @param pheno A matrix of phenotypes, individuals x phenotypes.
#' @param addcovar An optional matrix of additive covariates.
#' @param Xcovar An optional matrix with additional additive covariates used for
#' null hypothesis when scanning the X chromosome.
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have \code{names}
#' for individual identifiers.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#' @param ... Additional control parameters; see Details.
#'
#' @return A matrix of LOD scores, positions x phenotypes.
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals. The \code{genoprobs}
#' object should have an attribute \code{"is_x_chr"} that indicates
#' which of the chromosomes is the X chromosome, if any.
#'
#' The \code{...} argument can contain two additional control
#' parameters; suspended for simplicity (or confusion, depending on
#' your point of view). \code{"tol"} is used as a tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted);
#' default \code{1e-12}. \code{"intcovar_method"} indicates whether to
#' use a high-memory (but potentially faster) method or a low-memory
#' (and possibly slower) method, with values \code{"highmem"} or
#' \code{"lowmem"}; default \code{"lowmem"}.
#'
#' @references Haley, C. S. and Knott, S. A. (1992) A simple
#' regression method for mapping quantitative trait loci in line
#' crosses using flanking markers.  \emph{Heredity} \bold{69},
#' 315--324.
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
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, covar, Xcovar)
#'
#' @export
scan1 <-
    function(genoprobs, pheno, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, cores=1, ...)
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

    # force things to be matrices
    if(!is.matrix(pheno))
        pheno <- as.matrix(pheno)
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(Xcovar) && !is.matrix(Xcovar))
        Xcovar <- as.matrix(Xcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)
    # square-root of weights
    weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs[[1]], addcovar, Xcovar, intcovar,
                               weights, complete.cases=TRUE)
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
    if("cluster" %in% class(cores) && "SOCKcluster" %in% class(cores)) { # cluster already set
        cluster_ready <- TRUE
        n_cores <- length(cores)
        if(!quiet) message(" - Using ", n_cores, " cores.")
        quiet <- TRUE # no more messages
    } else {
        cluster_ready <- FALSE
        if(cores==0) cores <- parallel::detectCores() # if 0, detect cores
        if(cores > 1) {
            if(!quiet) message(" - Using ", cores, " cores.")
            quiet <- TRUE # no more messages

            if(Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
                n_cores <- cores
                cores <- parallel::makeCluster(cores)
                cluster_ready <- TRUE
                on.exit(parallel::stopCluster(cores))
            }
        }
        n_cores <- cores
    }

    # batches for analysis, to allow parallel analysis
    run_batches <- data.frame(chr=rep(seq(along=genoprobs), length(phe_batches)),
                              phe_batch=rep(seq(along=phe_batches), each=length(genoprobs)))
    run_indexes <- 1:(length(genoprobs)*length(phe_batches))

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
        ac <- addcovar; if(!is.null(ac)) ac <- ac[these2keep,,drop=FALSE]
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) ic <- ic[these2keep,,drop=FALSE]
        ph <- pheno[these2keep,phecol,drop=FALSE]
        wts <- weights[these2keep]

        # FIX_ME: calculating null RSS multiple times :(
        nullrss <- nullrss_clean(ph, ac, wts, tol)

        # if X chr, paste X covariates onto additive covariates
        if(is_x_chr[chr]) ac <- cbind(ac, Xc)

        # scan1 function taking clean data (with no missing values)
        rss <- scan1_clean(pr, ph, ac, ic, wts, tol, intcovar_method)

        # calculate LOD score
        nrow(ph)/2 * (log10(nullrss) - log10(rss))
    }

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- vapply(genoprobs, function(a) dim(a)[3], 1)
    totpos <- sum(npos_by_chr)
    pos_index <- split(1:totpos, rep(seq(along=genoprobs), npos_by_chr))

    # object to contain the LOD scores
    result <- matrix(nrow=totpos, ncol=ncol(pheno))

    if(cores<=1) { # no parallel processing
        for(i in run_indexes) {
            chr <- run_batches$chr[i]
            chrnam <- names(genoprobs)[chr]
            phebatch <- phe_batches[[run_batches$phe_batch[i]]]
            phecol <- phebatch$cols

            tmp_result <- by_group_func(i)
            if(!is.null(tmp_result))
                result[pos_index[[chr]], phecol] <- t(tmp_result)
        }
    }
    else {
        # calculations in parallel
        if(cluster_ready) # Windows doesn't suport mclapply
            list_result <- parallel::clusterApply(cores, run_indexes, by_group_func)
        else
            list_result <- parallel::mclapply(run_indexes, by_group_func, mc.cores=cores)

        # reorganize results
        for(i in run_indexes) {
            chr <- run_batches$chr[i]
            chrnam <- names(genoprobs)[chr]
            phebatch <- phe_batches[[run_batches$phe_batch[i]]]
            phecol <- phebatch$cols

            if(!is.null(list_result[[i]]))
                result[pos_index[[chr]], phecol] <- t(list_result[[i]])
        }
    }

    pos_names <- unlist(lapply(genoprobs, function(a) dimnames(a)[[3]]))
    dimnames(result) <- list(pos_names, colnames(pheno))

    # add some attributes with details on analysis
    attr(result, "addcovar") <- colnames4attr(addcovar)
    attr(result, "Xcovar") <- colnames4attr(Xcovar)
    attr(result, "intcovar") <- colnames4attr(intcovar)
    if(!is.null(weights))
        attr(result, "weights") <- TRUE

    result
}


# scan1 function taking nicely aligned data with no missing values
scan1_clean <-
    function(genoprobs, pheno, addcovar, intcovar,
             weights, tol, intcovar_method)
{
    n <- nrow(pheno)
    addcovar <- cbind(rep(1,n), addcovar) # add intercept

    if(is.null(intcovar)) { # no interactive covariates

        if(is.null(weights)) { # no weights
            return( scan_hk_onechr(genoprobs, pheno, addcovar, tol) )
        } else { # weights included
            return( scan_hk_onechr_weighted(genoprobs, pheno, addcovar, weights, tol) )
        }

    } else { # interactive covariates
        # high- and low-memory versions of functions
        if(intcovar_method=="highmem")
            scanf <- c(scan_hk_onechr_intcovar_highmem,
                       scan_hk_onechr_intcovar_weighted_highmem)
        else
            scanf <- c(scan_hk_onechr_intcovar_lowmem,
                       scan_hk_onechr_intcovar_weighted_lowmem)

        if(is.null(weights)) { # no weights
            return( scanf[[1]](genoprobs, pheno, addcovar, intcovar, tol) )
        } else { # weights included
            return( scanf[[2]](genoprobs, pheno, addcovar, intcovar, weights, tol) )
        }

    }
}

# calculate null RSS, with nicely aligned data with no missing values
nullrss_clean <-
    function(pheno, addcovar, weights, tol)
{
    n <- nrow(pheno)
    addcovar <- cbind(rep(1,n), addcovar) # add intercept

    if(is.null(weights)) { # no weights
        result <- calc_rss_linreg(addcovar, pheno, tol)
    } else { # weights included
        result <- calc_rss_linreg(addcovar*weights, pheno*weights, tol)
    }

    as.numeric(result)
}

# function to add column names as attribute
colnames4attr <-
    function(mat)
{
    if(is.null(mat)) return(mat)

    if(!is.matrix(mat)) mat <- as.matrix(mat)
    cn <- colnames(mat)

    if(is.null(cn)) return(character(ncol(mat)))

    cn
}
