#' Calculate QTL effects in scan along one chromosome
#'
#' Calculate QTL effects in scan along one chromosome with a
#' single-QTL model using Haley-Knott regression, with possible
#' allowance for covariates.
#'
#' @param genoprobs A 3-dimensional array of genotype probabilities
#' with dimension individuals x genotypes x positions.
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param addcovar An optional matrix of additive covariates.
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have \code{names}
#' for individual identifiers.
#' @param contrasts An optional matrix of genotype contrasts, size
#' genotypes x genotypes. For an intercross, you might use
#' \code{cbind(c(1,0,0), c(-0.5, 0, 0.5), c(-0.5, 1, 0.5))} to get
#' mean, additive effect, and dominance effect. The default is the
#' identity matrix.
#' @param se If TRUE, also calculate the standard errors.
#' @param tol Tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted)
#'
#' @return A matrix of estimated regression coefficients, of dimension
#' positions x number of effects. The number of effects is
#' \code{n_genotypes + n_addcovar + (n_genotypes-1)*n_intcovar}. The
#' map of positions at which the calculations were performed is
#' included as an attribute \code{"map"} (taken from the corresponding
#' attribute in the input \code{genoprobs}).
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals.
#'
#' @references Haley CS, Knott SA (1992) A simple
#' regression method for mapping quantitative trait loci in line
#' crosses using flanking markers.  Heredity 69:315--324.
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
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' coef <- scan1coef(probs[,7], pheno, covar)
#'
#' @export
scan1coef <-
    function(genoprobs, pheno, addcovar=NULL, intcovar=NULL, weights=NULL,
             contrasts=NULL, se=FALSE, tol=1e-12)
{
    stopifnot(tol > 0)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)
    if(!is.null(contrasts) && !is.matrix(contrasts))
        contrasts <- as.matrix(contrasts)
    # square-root of weights
    weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)

    # make sure pheno is a vector
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno) > 1)
            warning("Considering only the first phenotype.")
        rn <- rownames(pheno)
        pheno <- pheno[,1]
        names(pheno) <- rn
    }

    # genoprobs a list? then just take first chromosome
    if(is.list(genoprobs)) {
        if(length(genoprobs) > 1)
            warning("Considering only the first chromosome.")
        map <- attr(genoprobs, "map")[[1]]
        genoprobs <- genoprobs[[1]]
    } else {
        map <- attr(genoprobs, "map")
    }

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids(genoprobs, pheno, addcovar, intcovar,
                               weights, complete.cases=TRUE)
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
    if(!is.null(intcovar)) intcovar <- intcovar[ind2keep,,drop=FALSE]
    if(!is.null(weights)) weights <- weights[ind2keep]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # if weights, adjust phenotypes
    if(!is.null(weights)) pheno <- weights * pheno

    # weights have 0 dimension if missing
    if(is.null(weights)) weights <- numeric(0)

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    if(se) { # also calculate SEs

        if(is.null(addcovar))      # no covariates
            result <- scancoefSE_hk_nocovar(genoprobs, pheno, weights, tol)
        else if(is.null(intcovar)) # just addcovar
            result <- scancoefSE_hk_addcovar(genoprobs, pheno, addcovar, weights, tol)
        else                       # intcovar
            result <- scancoefSE_hk_intcovar(genoprobs, pheno, addcovar, intcovar,
                                             weights, tol)

        # move SEs to attribute
        SE <- t(result$SE) # transpose to positions x coefficients
        result <- result$coef
        attr(result, "SE") <- SE

    } else { # don't calculate SEs

        if(is.null(addcovar))      # no covariates
            result <- scancoef_hk_nocovar(genoprobs, pheno, weights, tol)
        else if(is.null(intcovar)) # just addcovar
            result <- scancoef_hk_addcovar(genoprobs, pheno, addcovar, weights, tol)
        else                       # intcovar
            result <- scancoef_hk_intcovar(genoprobs, pheno, addcovar, intcovar,
                                             weights, tol)
    }

    result <- t(result) # transpose to positions x coefficients

    # add names
    dimnames(result) <- list(dimnames(genoprobs)[[3]],
                             scan1coef_names(genoprobs, addcovar, intcovar))
    if(se) {
        dimnames(SE) <- dimnames(result)
        attr(result, "SE") <- SE
    }

    # add some attributes with details on analysis
    attr(result, "map") <- map
    attr(result, "sample_size") <- length(ind2keep)
    attr(result, "addcovar") <- colnames4attr(addcovar)
    attr(result, "intcovar") <- colnames4attr(intcovar)
    attr(result, "contrasts") <- contrasts
    if(!is.null(weights))
        attr(result, "weights") <- TRUE

    result
}


# genoprob x contrasts
genoprobs_by_contrasts <-
    function(genoprobs, contrasts)
{
    dg <- dim(genoprobs)
    dc <- dim(contrasts)
    if(dc[1] != dc[2] || dc[1] != dg[2])
        stop("contrasts should be a square matrix, ", dg[2], " x ", dg[2])

    # rearrange to put genotypes in last position
    dn <- dimnames(genoprobs)
    genoprobs <- aperm(genoprobs, c(1,3,2))
    dim(genoprobs) <- c(dg[1]*dg[3], dg[2])

    # multiply by contrasts
    genoprobs <- genoprobs %*% contrasts
    dim(genoprobs) <- dg[c(1,3,2)]

    genoprobs <- aperm(genoprobs, c(1,3,2))
    dn[[2]] <- colnames(contrasts)
    dimnames(genoprobs) <- dn

    genoprobs
}

# coefficient names
scan1coef_names <-
    function(genoprobs, addcovar, intcovar)
{
    qtl_names <-colnames(genoprobs)
    if(is.null(qtl_names))
        qtl_names <- paste0("qtleff", 1:ncol(genoprobs))

    if(is.null(addcovar)) { # no additive covariates
        return(qtl_names)
    }
    else { # some additive covariates
        add_names <- colnames(addcovar)
        if(is.null(add_names))
            add_names <- paste0("ac", 1:ncol(addcovar))

        if(is.null(intcovar)) { # no interactive covariates
            return(c(qtl_names, add_names))
        }

        int_names <- colnames(intcovar)
        if(is.null(intcovar))
            int_names <- paste0("ic", 1:ncol(intcovar))

        result <- c(qtl_names, add_names)
        for(ic in int_names)
            result <- c(result, paste(qtl_names[-1], ic, sep=":"))

        return(result)
    }
}
