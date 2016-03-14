#' Calculate BLUPs of QTL effects in scan along one chromosome
#'
#' Calculate BLUPs of QTL effects in scan along one chromosome, with a
#' single-QTL model treating the QTL effects as random, with possible
#' allowance for covariates and for a residual polygenic effect.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' \code{\link[qtl2geno]{calc_genoprob}}.
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param kinship Optional kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional matrix of additive covariates.
#' @param contrasts An optional matrix of genotype contrasts, size
#' genotypes x genotypes. For an intercross, you might use
#' \code{cbind(c(1,0,0), c(-0.5, 0, 0.5), c(-0.5, 1, 0.5))} to get
#' mean, additive effect, and dominance effect. The default is the
#' identity matrix.
#' @param se If TRUE, also calculate the standard errors.
#' @param reml If \code{reml=TRUE}, use
#' REML to estimate variance components; otherwise maximum likelihood.
#' @param preserve_intercept If TRUE, the BLUPs will have mean zero
#'     and there will be a separate column for the intercept. If FALSE
#'     (the default), we'll add the intercept to the BLUPs to give results
#'     that are comparable to \code{\link{scan1coef}}.
#' @param tol Tolerance value for convergence of linear mixed model fit.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A list containing the following
#' \itemize{
#' \item \code{coef} - A matrix of estimated regression coefficients, of dimension
#'     positions x number of effects. The number of effects is
#'     \code{n_genotypes + n_addcovar + (n_genotypes-1)*n_intcovar}.
#' \item \code{map} - A list containing the map positions at which the
#'     calculations were performed, taken from the input \code{genoprobs}.
#' \item \code{SE} - Present if \code{se=TRUE}: a matrix of estimated
#'     standard errors, of same dimension as \code{coef}.
#' \item \code{contrasts} - The input matrix of genotype coefficient
#'     contrasts that were used.
#' \item \code{addcovar} - Names of additive covariates that were used.
#' }
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals.
#'
#' If \code{kinship} is provided, the linear mixed model accounts for
#' a residual polygenic effect, with a the polygenic variance
#' estimated under the null hypothesis of no (major) QTL, and then
#' taken as fixed as known in the scan to estimate QTL effects.
#'
#' @references Haley CS, Knott SA (1992) A simple
#' regression method for mapping quantitative trait loci in line
#' crosses using flanking markers.  Heredity 69:315--324.
#'
#' Robinson GK (1991) That BLUP is a good thing: The estimation of
#' random effects. Statist Sci 6:15--32.
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
#' # convert to allele probabilities
#' aprobs <- genoprob_to_alleleprob(probs)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate BLUPs of coefficients for chromosome 7
#' blup <- scan1blup(aprobs[,7], pheno, addcovar=covar)
#'
#' # leave-one-chromosome-out kinship matrix for chr 7
#' kinship7 <- calc_kinship(probs, "loco")[[7]]
#'
#' # calculate BLUPs of coefficients for chromosome 7, adjusting for residual polygenic effect
#' blup_pg <- scan1blup(aprobs[,7], pheno, kinship7, addcovar=covar)
#'
#' @export
scan1blup <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL,
             contrasts=NULL, se=FALSE, reml=TRUE, preserve_intercept=FALSE,
             tol=1e-12, cores=1, quiet=TRUE)
{
    if(!is.null(kinship)) { # use LMM; see scan1_pg.R
        return(scan1blup_pg(genoprobs, pheno, kinship, addcovar,
                            contrasts, se, reml, preserve_intercept, tol, cores, quiet))
    }

    stopifnot(tol > 0)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(contrasts) && !is.matrix(contrasts))
        contrasts <- as.matrix(contrasts)

    # make sure pheno is a vector
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno) > 1)
            warning("Considering only the first phenotype.")
        rn <- rownames(pheno)
        pheno <- pheno[,1]
        names(pheno) <- rn
    }

    # genoprobs has more than one chromosome?
    if(length(genoprobs$probs) > 1)
        warning("Using only the first chromosome, ", names(genoprobs)[1])
    map <- genoprobs$map[[1]]
    genoprobs <- genoprobs$probs[[1]]

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, addcovar, complete.cases=TRUE)
    if(length(ind2keep)<=2) {
        if(length(ind2keep)==0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse=":"))
    }

    # omit individuals not in common
    genoprobs <- genoprobs[ind2keep,,,drop=FALSE]
    pheno <- pheno[ind2keep]
    if(!is.null(addcovar)) addcovar <- addcovar[ind2keep,,drop=FALSE]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)
    # add the intercept
    addcovar <- cbind(intercept=rep(1, length(ind2keep)),
                      addcovar)
    rownames(addcovar) <- ind2keep

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # determine batches
    n_pos <- dim(genoprobs)[[3]]
    batches <- batch_vec(1:n_pos, ceiling(n_pos/n_cores(cores)))

    by_group_func <- function(i)
        scanblup(genoprobs[,,batches[[i]],drop=FALSE], pheno, addcovar, se, reml, preserve_intercept, tol)

    # scan to get BLUPs and coefficient estimates
    if(n_cores(cores)==1) {
        result <- scanblup(genoprobs, pheno, addcovar, se, reml, preserve_intercept, tol)
        coef <- t(result$coef)
        if(se) SE <- t(result$SE)
        else SE <- NULL
    } else {
        message("Multi-core")
        result <- cluster_lapply(cores, seq(along=batches), by_group_func)
        SE <- coef <- matrix(nrow=n_pos,ncol=nrow(result[[1]]$coef))
        for(i in seq(along=result)) {
            coef[batches[[i]],] <- t(result[[i]]$coef)
            if(se) SE[batches[[i]],] <- t(result[[i]]$SE)
            else SE <- NULL
        }
    }

    # add names
    if(preserve_intercept)
        coefnames <- scan1coef_names(genoprobs, addcovar, NULL)
    else
        coefnames <- scan1coef_names(genoprobs, addcovar[,-1,drop=FALSE], NULL)
    dimnames(coef) <- list(dimnames(genoprobs)[[3]], coefnames)
    if(se) dimnames(SE) <- dimnames(coef)

    result <- list(coef = coef,
                   map = map,
                   sample_size = length(ind2keep),
                   addcovar = colnames4attr(addcovar),
                   contrasts = contrasts)
    result$SE <- SE # include only if not NULL

    class(result) <- c("scan1coef", "scan1", "matrix")
    result
}
