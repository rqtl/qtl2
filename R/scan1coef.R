#' Calculate QTL effects in scan along one chromosome
#'
#' Calculate QTL effects in scan along one chromosome with a
#' single-QTL model using Haley-Knott regression or a linear mixed
#' model (the latter to account for a residual polygenic effect), with
#' possible allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param kinship Optional kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional numeric matrix of additive covariates.
#' @param nullcovar An optional numeric matrix of additional additive
#' covariates that are used under the null hypothesis (of no QTL) but
#' not under the alternative (with a QTL). This is needed for the X
#' chromosome, where we might need sex as a additive covariate under
#' the null hypothesis, but we wouldn't want to include it under the
#' alternative as it would be collinear with the QTL effects. Only
#' used if `kinship` is provided but `hsq` is not, to get
#' estimate of residual heritability.
#' @param intcovar An optional numeric matrix of interactive covariates.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param contrasts An optional numeric matrix of genotype contrasts, size
#' genotypes x genotypes. For an intercross, you might use
#' `cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0))` to get
#' mean, additive effect, and dominance effect. The default is the
#' identity matrix.
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param zerosum If TRUE, force the genotype or allele coefficients
#'     sum to 0 by subtracting their mean and add another column with
#'     the mean. Ignored if `contrasts` is provided.
#' @param se If TRUE, also calculate the standard errors.
#' @param hsq (Optional) residual heritability; used only if
#' `kinship` provided.
#' @param reml If `kinship` provided: if `reml=TRUE`, use
#' REML; otherwise maximum likelihood.
#' @param ... Additional control parameters; see Details;
#'
#' @return A matrix of estimated regression coefficients, of dimension
#'     positions x number of effects. The number of effects is
#'     `n_genotypes + n_addcovar + (n_genotypes-1)*n_intcovar`.
#' May also contain the following attributes:
#' * `SE` - Present if `se=TRUE`: a matrix of estimated
#'   standard errors, of same dimension as `coef`.
#' * `sample_size` - Vector of sample sizes used for each
#'   phenotype
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals.
#'
#' If `kinship` is absent, Haley-Knott regression is performed.
#' If `kinship` is provided, a linear mixed model is used, with a
#' polygenic effect estimated under the null hypothesis of no (major)
#' QTL, and then taken as fixed as known in the genome scan.
#'
#' If `contrasts` is provided, the genotype probability matrix,
#' \eqn{P}, is post-multiplied by the contrasts matrix, \eqn{A}, prior
#' to fitting the model. So we use \eqn{P \cdot A}{P A} as the \eqn{X}
#' matrix in the model. One might view the rows of
#' \ifelse{html}{\out{<em>A</em><sup>-1</sup>}}{\eqn{A^{-1}}}
#' as the set of contrasts, as the estimated effects are the estimated
#' genotype effects pre-multiplied by
#' \ifelse{html}{\out{<em>A</em><sup>-1</sup>}}{\eqn{A^{-1}}}.
#'
#' The `...` argument can contain several additional control
#' parameters; suspended for simplicity (or confusion, depending on
#' your point of view). `tol` is used as a tolerance value for linear
#' regression by QR decomposition (in determining whether columns are
#' linearly dependent on others and should be omitted); default
#' `1e-12`. `maxit` is the maximum number of iteractions for
#' converence of the iterative algorithm used when `model=binary`.
#' `bintol` is used as a tolerance for converence for the iterative
#' algorithm used when `model=binary`. `eta_max` is the maximum value
#' for the "linear predictor" in the case `model="binary"` (a bit of a
#' technicality to avoid fitted values exactly at 0 or 1).
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
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#' \dontshow{probs[["7"]] <- probs[["7"]][,,1:5] # reduce to very small number}
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' coef <- scan1coef(probs[,7], pheno, addcovar=covar)
#'
#' # leave-one-chromosome-out kinship matrix for chr 7
#' kinship7 <- calc_kinship(probs, "loco")[[7]]
#'
#' # calculate coefficients for chromosome 7, adjusting for residual polygenic effect
#' coef_pg <- scan1coef(probs[,7], pheno, kinship7, addcovar=covar)
#'
#'
#' @export
scan1coef <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, nullcovar=NULL,
             intcovar=NULL, weights=NULL,
             contrasts=NULL, model=c("normal", "binary"), zerosum=TRUE,
             se=FALSE, hsq=NULL, reml=TRUE, ...)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    if(!is.null(kinship)) { # use LMM; see scan1_pg.R
        return(scan1coef_pg(genoprobs, pheno, kinship, addcovar, nullcovar,
                            intcovar, weights, contrasts, zerosum, se, hsq, reml, ...))
    }

    model <- match.arg(model)

    # deal with the dot args
    dotargs <- list("...")
    tol <- grab_dots(dotargs, "tol", 1e-12)
    if(!is_pos_number(tol)) stop("tol should be a single positive number")
    if(model=="binary") {
        bintol <- grab_dots(dotargs, "bintol", sqrt(tol)) # for model="binary"
        if(!is_pos_number(bintol)) stop("bintol should be a single positive number")
        eta_max <- grab_dots(dotargs, "eta_max", log(1-tol)-log(tol)) # for model="binary"
        if(!is_pos_number(eta_max)) stop("eta_max should be a single positive number")
        maxit <- grab_dots(dotargs, "maxit", 100) # for model="binary"
        if(!is_nonneg_number(maxit)) stop("maxit should be a single non-negative integer")

        check_extra_dots(dotargs, c("tol", "bintol", "eta_max", "maxit"))
    }
    else { # not binary trait
        check_extra_dots(dotargs, "tol")
    }

    # check that the objects have rownames
    check4names(pheno, addcovar, NULL, intcovar, nullcovar)

    # force things to be matrices
    if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.numeric(addcovar)) stop("addcovar is not numeric")
    }
    if(!is.null(nullcovar)) {
        if(!is.matrix(nullcovar)) nullcovar <- as.matrix(nullcovar)
        if(!is.numeric(nullcovar)) stop("nullcovar is not numeric")
    }
    if(!is.null(intcovar)) {
        if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)
        if(!is.numeric(intcovar)) stop("intcovar is not numeric")
    }
    if(!is.null(contrasts)) {
        if(!is.matrix(contrasts)) contrasts <- as.matrix(contrasts)
        if(!is.numeric(contrasts)) stop("contrasts is not numeric")
    }

    # make sure pheno is a vector
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno) > 1)
            warning("Considering only the first phenotype.")
        rn <- rownames(pheno)
        pheno <- pheno[,1]
        names(pheno) <- rn
        if(!is.numeric(pheno)) stop("pheno is not numeric")
    }

    # genoprobs has more than one chromosome?
    if(length(genoprobs) > 1)
        warning("Using only the first chromosome, ", names(genoprobs)[1])
    chrid <- names(genoprobs)[1]
    genoprobs <- genoprobs[[1]]

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

    # normal or binary model?
    if(model=="binary") {
        if(!is.null(kinship))
            stop("Can't yet account for kinship with model = \"binary\"")
        if(any(!is.na(pheno) & (pheno < 0 | pheno > 1)))
            stop('with model="binary", pheno should be in [0,1]')
    }
    else {
        # square-root of weights (only if model="normal")
        weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)

        # if weights, adjust phenotypes
        if(!is.null(weights)) pheno <- weights * pheno
    }

    # weights have 0 dimension if missing
    if(is.null(weights)) weights <- numeric(0)

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs_by_contrasts(genoprobs, contrasts)

    if(se) { # also calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            if(model=="normal")
                result <- scancoefSE_hk_addcovar(genoprobs, pheno, addcovar, weights, tol)
            else # binary trait
                result <- scancoefSE_binary_addcovar(genoprobs, pheno, addcovar, weights, maxit, bintol, tol, eta_max)
        }
        else {                  # intcovar
            if(model=="normal")
                result <- scancoefSE_hk_intcovar(genoprobs, pheno, addcovar, intcovar,
                                                 weights, tol)
            else
                result <- scancoefSE_binary_intcovar(genoprobs, pheno, addcovar, intcovar,
                                                 weights, maxit, bintol, tol, eta_max)
        }

        SE <- t(result$SE) # transpose to positions x coefficients
        result <- result$coef
    } else { # don't calculate SEs

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            if(model=="normal")
                result <- scancoef_hk_addcovar(genoprobs, pheno, addcovar, weights, tol)
            else
                result <- scancoef_binary_addcovar(genoprobs, pheno, addcovar, weights, maxit, bintol, tol, eta_max)
        }
        else {                  # intcovar
            if(model=="normal")
                result <- scancoef_hk_intcovar(genoprobs, pheno, addcovar, intcovar,
                                               weights, tol)
            else
                result <- scancoef_binary_intcovar(genoprobs, pheno, addcovar, intcovar,
                                                   weights, maxit, bintol, tol, eta_max)
        }
        SE <- NULL
    }

    result <- t(result) # transpose to positions x coefficients

    # add names
    dimnames(result) <- list(dimnames(genoprobs)[[3]],
                             scan1coef_names(genoprobs, addcovar, intcovar))
    if(se) dimnames(SE) <- dimnames(result)

    if(zerosum && is.null(contrasts)) { # force QTL effects to sum to 0
        ng <- dim(genoprobs)[2]
        whcol <- seq_len(ng)
        mu <- rowMeans(result[,whcol,drop=FALSE], na.rm=TRUE)
        result <- cbind(result, mean=mu)
        result[,whcol] <- result[,whcol] - mu

        if(se) {
            SE <- cbind(SE, mean=sqrt(rowMeans(SE[,whcol,drop=FALSE]^2)))
        }

    }

    attr(result, "sample_size") <- length(ind2keep)
    attr(result, "SE") <- SE # include only if not NULL

    class(result) <- c("scan1coef", "scan1", "matrix")
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
        qtl_names <- paste0("qtleff", seq_len(ncol(genoprobs)))

    if(is.null(addcovar) || ncol(addcovar)==0) { # no additive covariates
        return(qtl_names)
    }
    else { # some additive covariates
        add_names <- colnames(addcovar)
        if(is.null(add_names) || all(add_names==""))
            add_names <- paste0("ac", seq_len(ncol(addcovar)))
        else if(all(add_names[-1] == "")) # all but first is empty
            add_names[-1] <- paste0("ac", seq_len(ncol(addcovar)-1))

        if(is.null(intcovar)) { # no interactive covariates
            return(c(qtl_names, add_names))
        }

        int_names <- colnames(intcovar)
        if(is.null(int_names) || all(int_names==""))
            int_names <- paste0("ic", seq_len(ncol(intcovar)))
        else if(all(int_names[-1] == "")) # all but first is empty
            int_names[-1] <- paste0("ic", seq_len(ncol(intcovar)-1))

        result <- c(qtl_names, add_names)
        for(ic in int_names)
            result <- c(result, paste(qtl_names[-1], ic, sep=":"))

        return(result)
    }
}
