#' Permutation test for genome scan with a single-QTL model
#'
#' Permutation test for a enome scan with a single-QTL model by
#' Haley-Knott regression or a linear mixed model, with possible
#' allowance for covariates.
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
#' @param n_perm Number of permutation replicates.
#' @param perm_Xsp If TRUE, do separate permutations for the autosomes
#' and the X chromosome.
#' @param perm_strata Vector of strata, for a stratified permutation
#' test. Should be named in the same way as the rows of
#' \code{pheno}. The unique values define the strata.
#' @param chr_lengths Lengths of the chromosomes; needed only if
#' \code{perm_Xsp=TRUE}. See \code{\link{chr_lengths}}.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#' @param ... Additional control parameters; see Details.
#'
#' @return If \code{perm_Xsp=FALSE}, the result is matrix of
#' genome-wide maximum LOD scores, permutation replicates x
#' phenotypes. If \code{perm_Xsp=TRUE}, the result is a list of
#' two matrices, one for the autosomes and one for the X
#' chromosome.
#'
#' @details
#' If \code{kinship} is not provided, so that analysis proceeds by
#' Haley-Knott regression, we permute the rows of the phenotype data;
#' the same permutations are also applied to the rows of the
#' covariates (\code{addcovar}, \code{Xcovar}, and \code{intcovar})
#' are permuted.
#'
#' If \code{kinship} is provided, we instead permute the rows of the
#' genotype data and fit an LMM with the same residual heritability
#' (estimated under the null hypothesis of no QTL).
#'
#' If \code{Xcovar} is provided and \code{perm_strata=NULL}, we do a
#' stratified permutation test with the strata defined by the rows of
#' \code{Xcovar}. If a simple permutation test is desired, provide
#' \code{perm_strata} that is a vector containing a single repeated
#' value.
#'
#' The \code{...} argument can contain three additional control
#' parameters; suspended for simplicity (or confusion, depending on
#' your point of view). \code{tol} is used as a tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted);
#' default \code{1e-12}. \code{intcovar_method} indicates whether to
#' use a high-memory (but potentially faster) method or a low-memory
#' (and possibly slower) method, with values \code{"highmem"} or
#' \code{"lowmem"}; default \code{"lowmem"}.  Finally, \code{max_batch}
#' indicates the maximum number of phenotypes to run together; default
#' is 1000.
#'
#' @references Churchill GA, Doerge RW (1994) Empirical threshold
#' values for quantitative trait mapping. Genetics 138:963--971.
#'
#' Manichaikul A, Palmer AA, Sen S, Broman KW (2007) Significance
#' thresholds for quantitative trait locus mapping under selective
#' genotyping. Genetics 177:1963--1966.
#'
#' Haley CS, Knott SA (1992) A simple regression method for mapping
#' quantitative trait loci in line crosses using flanking markers.
#' Heredity 69:315--324.
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
#' # strata for permutations
#' perm_strata <- mat2strata(Xcovar)
#'
#' # permutations with genome scan
#' \dontrun{
#' operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
#'                    n_perm=1000, perm_Xsp=TRUE, perm_strata=perm_strata,
#'                    chr_lengths=chr_lengths(iron$gmap))}
#' \dontshow{operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)}
#'
#' # leave-one-chromosome-out kinship matrices
#' kinship <- calc_kinship(probs, "loco")
#'
#' # genome scan with a linear mixed model
#' \dontrun{
#' operm_lmm <- scan1perm(probs, pheno, kinship, covar, Xcovar, n_perm=1000,
#'                        perm_Xsp=TRUE, perm_strata=perm_strata)}
#' \dontshow{operm_lmm <- scan1perm(probs, pheno, kinship, covar, Xcovar, n_perm=3, perm_strata=perm_strata)}
#'
#' @seealso \code{\link{scan1}}, \code{\link{chr_lengths}}, \code{\link{mat2strata}}
#' @export
scan1perm <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, n_perm=1,
             perm_Xsp=FALSE, perm_strata=NULL, chr_lengths=NULL,
             cores=1, ...)
{
    # grab tol from dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    stopifnot(tol > 0)

    if(n_perm <= 0) stop("n_perm should be > 0")

    # if Xcovar provided, the default is to do a stratified permutation test
    if(is.null(perm_strata) && !is.null(Xcovar))
        perm_strata <- mat2strata(Xcovar)

    if(perm_Xsp) { # autosome/X chr-specific permutations
        if(is.null(chr_lengths))
            stop("Need to provide chr_lengths when perm_Xsp=TRUE")

        is_x_chr <- attr(genoprobs, "is_x_chr")

        # collapse to A/X
        chr_lengths <- collapse_chr_lengths_to_AX(chr_lengths, is_x_chr)
        if(any(chr_lengths < 1e-12))
            stop("Autosome or X chromosome of length 0; skip perm_Xsp")

        n_permX <- ceiling(n_perm * chr_lengths[["A"]] / chr_lengths[["X"]])

        A <- scan1perm(genoprobs=genoprobs[,!is_x_chr], pheno=pheno,
                       kinship=subset_kinship(kinship, chr=!is_x_chr),
                       addcovar=addcovar, Xcovar=NULL, intcovar=intcovar, weights=weights,
                       reml=reml, n_perm=n_perm, perm_Xsp=FALSE, perm_strata=perm_strata,
                       chr_lengths=NULL, cores=cores, ...)
        X <- scan1perm(genoprobs=genoprobs[,is_x_chr], pheno=pheno,
                       kinship=subset_kinship(kinship, chr=is_x_chr),
                       addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar, weights=weights,
                       reml=reml, n_perm=n_permX, perm_Xsp=FALSE, perm_strata=perm_strata,
                       chr_lengths=NULL, cores=cores, ...)
        result <- list(A=A, X=X)
        attr(result, "chr_lengths") <- chr_lengths

        class(result$A) <- class(result$X) <- "matrix"
        class(result) <- c("scan1perm", "list")
        return(result)
    }

    # force things to be matrices
    if(!is.matrix(pheno))
        pheno <- as.matrix(pheno)
    if(is.null(colnames(pheno))) # force column names
        colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(Xcovar) && !is.matrix(Xcovar))
        Xcovar <- as.matrix(Xcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)

    # check that kinship matrices are square with same IDs
    kinshipIDs <- check_kinship(kinship, length(genoprobs))

    # multiply kinship matrix by 2; rest is using 2*kinship
    # see Almasy & Blangero (1998) https://doi.org/10.1086/301844
    kinship <- double_kinship(kinship)

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
                               kinshipIDs, weights, perm_strata, complete.cases=TRUE)
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

    if(!is.null(kinship)) { # fit linear mixed model
        return(scan1perm_pg(genoprobs=genoprobs, pheno=pheno, kinship=kinship,
                            addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar,
                            reml=reml, n_perm=n_perm, perm_strata=perm_strata,
                            cores=cores, ind2keep=ind2keep, ...))
    }

    if(is.null(addcovar) && is.null(Xcovar) &&
       is.null(intcovar) && is.null(weights)
       && sum(is.na(pheno[ind2keep,]))==0) # no covariates, no weights, no missing phenotypes
        result <- scan1perm_nocovar(genoprobs=genoprobs,
                                    pheno=pheno,
                                    n_perm=n_perm,
                                    perm_strata=perm_strata,
                                    cores=cores,
                                    ind2keep=ind2keep,
                                    ...)
    else
        result <- scan1perm_covar(genoprobs=genoprobs,
                                  pheno=pheno,
                                  addcovar=addcovar,
                                  Xcovar=Xcovar,
                                  intcovar=intcovar,
                                  weights=weights,
                                  n_perm=n_perm,
                                  perm_strata=perm_strata,
                                  cores=cores,
                                  ind2keep=ind2keep,
                                  ...)


    result
}


# simplest version: no covariates, no weights, no missing phenotypes
scan1perm_nocovar <-
    function(genoprobs, pheno, n_perm=1, perm_strata=NULL,
             cores=1, ind2keep, ...)
{
    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    stopifnot(tol > 0)
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch"))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # max batch size
    max_batch <- grab_dots(dotargs, "max_batch",
                           min(1000, ceiling(n_perm*length(genoprobs)*ncol(pheno)/n_cores(cores))))

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

    nullrss <- apply(pheno, 2, function(ph) nullrss_clean(as.matrix(ph), NULL, NULL, add_intercept=TRUE, tol))

    # the function that does the work
    by_group_func <- function(i) {
        # deal with batch information, including individuals to drop due to missing phenotypes
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        permbatch <- perm_batches[[run_batches$phe_batch[i]]]

        # subset the genotype probabilities: drop cols with all 0s, plus the first column
        Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
        if(length(Xcol2drop) > 0) {
            pr <- genoprobs[[chr]][ind2keep,-Xcol2drop,,drop=FALSE]
            pr <- pr[,-1,,drop=FALSE]
        }
        else
            pr <- genoprobs[[chr]][ind2keep,-1,,drop=FALSE]

        ph <- pheno[,phebatch,drop=FALSE]
        for(col in seq_len(ncol(ph))) # permute columns
            ph[,col] <- ph[perms[,permbatch[col]] , col]

        # scan1 function taking clean data (with no missing values)
        rss <- scan1_clean(pr, ph, NULL, NULL, NULL, add_intercept=TRUE, tol, "lowmem")

        # calculate LOD score
        lod <- nind/2 * (log10(nullrss[phebatch]) - log10(rss))

        # return column maxima
        apply(lod, 1, max, na.rm=TRUE)
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

# permutations with covariates and/or different batches of phenotypes
scan1perm_covar <-
    function(genoprobs, pheno, addcovar=NULL, Xcovar=NULL, intcovar=NULL,
             weights=weights, n_perm=1, perm_strata=NULL, cores=1,
             ind2keep, ...)
{
    # deal with the dot args
    dotargs <- list(...)
    tol <- grab_dots(dotargs, "tol", 1e-12)
    stopifnot(tol > 0)
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch"))

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # max batch size
    max_batch <- grab_dots(dotargs, "max_batch",
                           min(1000, ceiling(n_perm*length(genoprobs)*ncol(pheno)/n_cores(cores))))

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
        pr <- pr[,-1,,drop=FALSE]

        # subset the rest
        ac <- addcovar; if(!is.null(ac)) ac <- ac[these2keep,,drop=FALSE]
        Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
        ic <- intcovar; if(!is.null(ic)) ic <- ic[these2keep,,drop=FALSE]
        wts <- weights; if(!is.null(wts)) wts <- wts[these2keep]
        ph <- pheno[these2keep, phecol, drop=FALSE]

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

        # return column maxima
        apply(lod, 1, max, na.rm=TRUE)
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
