#' Genome scan with a single-QTL model
#'
#' Genome scan with a single-QTL model by Haley-Knott regression or a
#' linear mixed model, with possible allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' \code{\link[qtl2geno]{calc_genoprob}}.
#' @param pheno A matrix of phenotypes, individuals x phenotypes.
#' @param kinship Optional kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional matrix of additive covariates.
#' @param Xcovar An optional matrix with additional additive covariates used for
#' null hypothesis when scanning the X chromosome.
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have \code{names}
#' for individual identifiers. Ignored if \code{kinship} is provided.
#' @param reml If \code{kinship} provided: if \code{reml=TRUE}, use
#' REML; otherwise maximum likelihood.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#' @param ... Additional control parameters; see Details.
#'
#' @return A matrix of LOD scores, positions x phenotypes.
#' Also contains one or more of the following attributes:
#' \itemize{
#' \item \code{sample_size} - Vector of sample sizes used for each
#'     phenotype
#' \item \code{hsq} - Included if \code{kinship} provided: A matrix of
#'     estimated heritabilities under the null hypothesis of no
#'     QTL. Columns are the phenotypes. If the \code{"loco"} method was
#'     used with \code{\link[qtl2geno]{calc_kinship}} to calculate a list
#'     of kinship matrices, one per chromosome, the rows of \code{hsq}
#'     will be the heritabilities for the different chromosomes (well,
#'     leaving out each one). If \code{Xcovar} was not NULL, there will at
#'     least be an autosome and X chromosome row.
#' }
#'
#' @details
#' We first fit the model \eqn{y = X \beta + \epsilon}{y = Xb + e}
#' where \eqn{X} is a matrix of covariates (or just an intercept) and
#' \eqn{\epsilon}{e} is multivariate normal with mean 0 and covariance
#' matrix \eqn{\sigma^2 [h^2 (2 K) + I]}{sigmasq*[hsq*2*K+I]} where
#' \eqn{K} is the kinship matrix and \eqn{I} is the identity matrix.
#'
#' We then take \eqn{h^2}{hsq} as fixed and then scan the genome, at
#' each genomic position fitting the model \eqn{y = P \alpha + X \beta
#' + \epsilon}{y = Xb + e} where \eqn{P} is a matrix of genotype
#' probabilities for the current position and again \eqn{X} is a
#' matrix of covariates \eqn{\epsilon}{e} is multivariate normal with
#' mean 0 and covariance matrix \eqn{\sigma^2 [h^2 (2 K) +
#' I]}{sigmasq*[hsq*2*K+I]}, taking \eqn{h^2}{hsq} to be known.
#'
#' For each of the inputs, the row names are used as
#' individual identifiers, to align individuals. The \code{genoprobs}
#' object should have a component \code{"is_x_chr"} that indicates
#' which of the chromosomes is the X chromosome, if any.
#'
#' The \code{...} argument can contain three additional control
#' parameters; suspended for simplicity (or confusion, depending on
#' your point of view). \code{tol} is used as a tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted);
#' default \code{1e-12}. \code{intcovar_method} indicates whether to
#' use a high-memory (but potentially faster) method or a low-memory
#' (and possibly slower) method, with values \code{"highmem"} or
#' \code{"lowmem"}; default \code{"lowmem"}. Finally, \code{max_batch}
#' indicates the maximum number of phenotypes to run together; default
#' is unlimited.
#'
#' If \code{kinship} is absent, Haley-Knott regression is performed.
#' If \code{kinship} is provided, a linear mixed model is used, with a
#' polygenic effect estimated under the null hypothesis of no (major)
#' QTL, and then taken as fixed as known in the genome scan.
#'
#' If \code{kinship} is a single matrix, then the \code{hsq}
#' in the results is a vector of heritabilities (one value for each phenotype). If
#' \code{kinship} is a list (one matrix per chromosome), then
#' \code{hsq} is a matrix, chromosomes x phenotypes.
#'
#' @references Haley CS, Knott SA (1992) A simple
#' regression method for mapping quantitative trait loci in line
#' crosses using flanking markers.  Heredity 69:315--324.
#'
#' Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, Eskin
#' E (2008) Efficient control of population structure in model
#' organism association mapping. Genetics 178:1709--1723.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
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
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # leave-one-chromosome-out kinship matrices
#' kinship <- calc_kinship(probs, "loco")
#'
#' # genome scan with a linear mixed model
#' out_lmm <- scan1(probs, pheno, kinship, covar, Xcovar)
#'
#' @seealso \code{\link{scan1perm}}
#' @export
scan1 <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, cores=1, ...)
{
    if(!is.null(kinship)) { # fit linear mixed model
        return(scan1_pg(genoprobs, pheno, kinship, addcovar, Xcovar, intcovar,
                        reml, cores, ...))
    }

    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    stopifnot(tol > 0)
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch"))

    # force things to be matrices
    if(!is.matrix(pheno))
        pheno <- as.matrix(pheno)
    if(is.null(colnames(pheno))) # force column names
        colnames(pheno) <- paste0("pheno", 1:ncol(pheno))
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
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
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
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
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

        # if X chr, paste X covariates onto additive covariates
        # (only for the null)
        if(is_x_chr[chr]) ac0 <- drop_depcols(cbind(ac, Xc), add_intercept=FALSE, tol)
        else ac0 <- ac

        # FIX_ME: calculating null RSS multiple times :(
        nullrss <- nullrss_clean(ph, ac0, wts, add_intercept=TRUE, tol)

        # scan1 function taking clean data (with no missing values)
        rss <- scan1_clean(pr, ph, ac, ic, wts, add_intercept=TRUE, tol, intcovar_method)

        # calculate LOD score
        lod <- nrow(ph)/2 * (log10(nullrss) - log10(rss))
        list(lod=lod, n=nrow(ph)) # return LOD & number of individuals used
    }

    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- vapply(genoprobs, function(a) dim(a)[3], 1)
    totpos <- sum(npos_by_chr)
    pos_index <- split(1:totpos, rep(seq(along=genoprobs), npos_by_chr))

    # object to contain the LOD scores; also attr to contain sample size
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)

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

    pos_names <- unlist(lapply(genoprobs, function(a) dimnames(a)[[3]]))
    names(pos_names) <- NULL # this is just annoying
    dimnames(result) <- list(pos_names, colnames(pheno))

    # add some attributes with details on analysis
    attr(result, "sample_size") <- n

    class(result) <- c("scan1", "matrix")
    result
}


# scan1 function taking nicely aligned data with no missing values
#
# Here genoprobs is a plain 3d array
scan1_clean <-
    function(genoprobs, pheno, addcovar, intcovar,
             weights, add_intercept=TRUE, tol, intcovar_method)
{
    n <- nrow(pheno)
    if(add_intercept)
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
    function(pheno, addcovar, weights, add_intercept=TRUE, tol)
{
    n <- nrow(pheno)
    if(add_intercept)
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
    if(is.null(mat)) return(NULL)

    if(!is.matrix(mat)) mat <- as.matrix(mat)
    if(ncol(mat)==0) return(NULL)
    cn <- colnames(mat)

    if(is.null(cn)) cn <- rep("", ncol(mat))

    if(any(cn=="")) cn[cn==""] <- paste0("unnamed", 1:sum(cn==""))

    cn
}
