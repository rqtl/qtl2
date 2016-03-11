#' Calculate QTL effects in scan along one chromosome, adjusting for polygenes with an LMM
#'
#' Calculate QTL effects in scan along one chromosome with a
#' single-QTL model using a linear mixed model, with possible
#' allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' \code{\link[qtl2geno]{calc_genoprob}}.
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param kinship A kinship matrix. Can be eigen decomposition of
#' kinship matrix (as calculated with \code{\link{decomp_kinship}}),
#' but it would need to be for the exact subset of individuals
#' (following omitting those with missing genotype probabilities,
#' phenotypes, or covariates).
#' @param addcovar An optional matrix of additive covariates.
#' @param intcovar An optional matrix of interactive covariates.
#' @param contrasts An optional matrix of genotype contrasts, size
#' genotypes x genotypes. For an intercross, you might use
#' \code{cbind(c(1,0,0), c(-0.5, 0, 0.5), c(-0.5, 1, 0.5))} to get
#' mean, additive effect, and dominance effect. The default is the
#' identity matrix.
#' @param se If TRUE, also calculate the standard errors.
#' @param hsq (Optional) residual heritability
#' @param reml If true and \code{hsq} is not provided, use REML to estimate \code{hsq}.
#' @param tol Tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted)
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
#' \item \code{intcovar} - Names of interactive covariates that were used.
#' \item \code{sample_size} - Vector of sample sizes used for each
#'     phenotype
#' }
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' \dontshow{iron <- iron[,c(7, 8)]}
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # leave-one-chromosome-out kinship matrices
#' kinship <- calc_kinship(probs, "loco")
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' coef <- scan1coef_lmm(probs[,"7"], pheno, kinship[["7"]], addcovar=covar)
#'
#' @export
scan1coef_lmm <-
    function(genoprobs, pheno, kinship,
             addcovar=NULL, intcovar=NULL,
             contrasts=NULL, se=FALSE,
             hsq=NULL, reml=TRUE, tol=1e-12)
{
    stopifnot(tol > 0)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)
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

    # check that kinship matrices are square with same IDs
    if(!is.null(attr(kinship, "eigen_decomp"))) { # already did decomposition
        kinshipIDs <- rownames(kinship$vectors)
        did_decomp <- TRUE
    } else {
        kinshipIDs <- check_kinship(kinship, 1)
        did_decomp <- FALSE
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, kinshipIDs,
                               addcovar, intcovar, complete.cases=TRUE)

    if(length(ind2keep)<=2) {
        if(length(ind2keep)==0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse=":"))
    }

    if(did_decomp) { # if did decomposition already, make sure it was with exactly
        if(length(kinshipIDs) != length(ind2keep) ||
           any(sort(kinshipIDs) != sort(ind2keep)))
            stop("Decomposed kinship matrix was with different individuals")
        else
            ind2keep <- kinshipIDs # force them in same order
    }

    # omit individuals not in common
    genoprobs <- genoprobs[ind2keep,,,drop=FALSE]
    pheno <- pheno[ind2keep]
    if(!is.null(addcovar)) addcovar <- addcovar[ind2keep,,drop=FALSE]

    if(!is.null(intcovar)) intcovar <- intcovar[ind2keep,,drop=FALSE]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # eigen decomposition of kinship matrix
    if(!did_decomp)
        kinship <- decomp_kinship(kinship[ind2keep, ind2keep])

    # estimate hsq if necessary
    if(is.null(hsq)) {
        nullresult <- calc_hsq_clean(kinship, as.matrix(pheno), addcovar, NULL, FALSE,
                                     reml, cores=1, check_boundary=TRUE, tol)
        hsq <- nullresult$hsq
    }

    # eigen-vectors and weights
    eigenvec <- kinship$vectors
    weights <- 1/sqrt(hsq*kinship$values + (1-hsq))

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    if(se) { # also calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            result <- scancoefSE_lmm_addcovar(genoprobs, pheno, addcovar, eigenvec, weights, tol)
        }
        else {                  # intcovar
            result <- scancoefSE_lmm_intcovar(genoprobs, pheno, addcovar, intcovar,
                                              eigenvec, weights, tol)
        }

        SE <- t(result$SE) # transpose to positions x coefficients
        result <- result$coef
    } else { # don't calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            result <- scancoef_lmm_addcovar(genoprobs, pheno, addcovar, eigenvec, weights, tol)
        }
        else {                  # intcovar
            result <- scancoef_lmm_intcovar(genoprobs, pheno, addcovar, intcovar,
                                            eigenvec, weights, tol)
        }
        SE <- NULL
    }

    result <- t(result) # transpose to positions x coefficients

    # add names
    dimnames(result) <- list(dimnames(genoprobs)[[3]],
                             scan1coef_names(genoprobs, addcovar, intcovar))
    if(se) dimnames(SE) <- dimnames(result)

    # add some attributes with details on analysis
    result <- list(coef = result,
                   map = map,
                   sample_size = length(ind2keep),
                   addcovar = colnames4attr(addcovar),
                   intcovar = colnames4attr(intcovar),
                   contrasts = contrasts)
    result$SE <- SE # include only if not NULL

    class(result) <- c("scan1coef", "scan1", "matrix")
    result
}
