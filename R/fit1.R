#' Fit single-QTL model at a single position
#'
#' Fit a single-QTL model at a single putative QTL position and get detailed results
#' about estimated coefficients and individuals contributions to the LOD score.
#'
#' @param genoprobs A matrix of genotype probabilities, individuals x genotypes
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param kinship Optional kinship matrix.
#' @param addcovar An optional numeric matrix of additive covariates.
#' @param nullcovar An optional numeric matrix of additional additive
#' covariates that are used under the null hypothesis (of no QTL)
#' but not under the alternative (with a QTL). This is needed for
#' the X chromosome, where we might need sex as a additive
#' covariate under the null hypothesis, but we wouldn't want to
#' include it under the alternative as it would be collinear with
#' the QTL effects.
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
#'     If `model="binary"`, the phenotypes must have values in
#'     \eqn{[0, 1]}.
#' @param zerosum If TRUE, force the genotype or allele coefficients
#'     sum to 0 by subtracting their mean and add another column with
#'     the mean. Ignored if `contrasts` is provided.
#' @param se If TRUE, calculate the standard errors.
#' @param hsq (Optional) residual heritability; used only if
#' `kinship` provided.
#' @param reml If `kinship` provided: if `reml=TRUE`, use
#' REML; otherwise maximum likelihood.
#' @param blup If TRUE, fit a model with QTL effects being random, as in [scan1blup()].
#' @param ... Additional control parameters; see Details;
#'
#' @return A list containing
#' * `coef` - Vector of estimated coefficients.
#' * `SE` - Vector of estimated standard errors (included if `se=TRUE`).
#' * `lod` - The overall lod score.
#' * `ind_lod` - Vector of individual contributions to the LOD score.
#' * `fitted`  - Fitted values.
#' If `blup==TRUE`, only `coef` and `SE` are included at present.
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
#' `1e-12`. `maxit` is the maximum number of iterations for
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
#' \dontshow{iron <- iron[,"7"] # subset to chr 7}
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=5)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # scan chromosome 7 to find peak
#' out <- scan1(probs[,"7"], pheno, addcovar=covar)
#'
#' # find peak position
#' max_pos <- max(out, map)
#'
#' # genoprobs at max position
#' pr_max <- pull_genoprobpos(probs, map, max_pos$chr, max_pos$pos)
#'
#' # fit QTL model just at that position
#' out_fit1 <- fit1(pr_max, pheno, addcovar=covar)
#'
#' @seealso [pull_genoprobpos()], [find_marker()]
#'
#' @importFrom stats setNames
#' @export
fit1 <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, nullcovar=NULL,
             intcovar=NULL, weights=NULL,
             contrasts=NULL, model=c("normal", "binary"),
             zerosum=TRUE, se=TRUE, hsq=NULL, reml=TRUE, blup=FALSE, ...)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    model <- match.arg(model)

    if(blup) {
        if(!is.null(intcovar)) warning("If blup==TRUE, intcovar ignored")
        if(model != "normal") warning('If blup==TRUE, model taken to be "normal"')
        if(!is.null(weights)) warning("If blup==TRUE, weights ignored")
        if(!is.null(hsq)) warning("If blup==TRUE, hsq ignored")

        rn <- rownames(genoprobs)
        cn <- colnames(genoprobs)
        genoprobs <- list("1" = array(genoprobs, dim=c(dim(genoprobs), 1)))
        dimnames(genoprobs[[1]]) <- list(rn, cn, "marker")

        blup_result <- scan1blup(genoprobs, pheno, kinship=kinship, addcovar=addcovar,
                                 nullcovar=nullcovar, contrasts=contrasts, se=se, reml=reml,
                                 tol=grab_dots(list("..."), "tol", 1e-12))
        result <- list(coef=blup_result[1,])
        if(se) result$SE <- attr(blup_result, "SE")[1,]
        return(result)
    }

    if(!is.null(kinship)) { # use LMM; see fit1_pg.R
        return(fit1_pg(genoprobs, pheno, kinship, addcovar, nullcovar,
                       intcovar, weights, contrasts, zerosum, se, hsq, reml, ...))
    }

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

    if(model=="binary") {
        if(!is.null(kinship))
            stop("Can't yet account for kinship with model = \"binary\"")
        if(any(!is.na(pheno) & (pheno < 0 | pheno > 1)))
            stop('with model="binary", pheno should be in [0,1]')
    }
    else {
        # square-root of weights
        weights <- sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)
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

    # genoprobs is a matrix?
    if(!is.matrix(genoprobs))
        stop("genoprobs should be a matrix, individuals x genotypes")

    # make sure contrasts is square n_genotypes x n_genotypes
    if(!is.null(contrasts)) {
        ng <- ncol(genoprobs)
        if(ncol(contrasts) != ng || nrow(contrasts) != ng)
            stop("contrasts should be a square matrix, ", ng, " x ", ng)
    }

    # find individuals in common across all arguments
    ind2keep <- get_common_ids(genoprobs, pheno, addcovar, nullcovar, intcovar,
                               weights, complete.cases=TRUE)
    if(length(ind2keep)<=2) {
        if(length(ind2keep)==0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse=":"))
    }

    # omit individuals not in common
    genoprobs <- genoprobs[ind2keep,,drop=FALSE]
    pheno <- pheno[ind2keep]
    if(!is.null(addcovar)) addcovar <- addcovar[ind2keep,,drop=FALSE]
    if(!is.null(nullcovar)) nullcovar <- nullcovar[ind2keep,,drop=FALSE]
    if(!is.null(intcovar)) intcovar <- intcovar[ind2keep,,drop=FALSE]
    if(!is.null(weights)) weights <- weights[ind2keep]

    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)

    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)

    # multiply genoprobs by contrasts
    if(!is.null(contrasts))
        genoprobs <- genoprobs %*% contrasts

    X0 <- drop_depcols(cbind(rep(1, length(pheno)), addcovar, nullcovar), FALSE, tol)
    if(model=="normal") {
        # if weights, adjust phenotypes
        if(!is.null(weights)) pheno <- weights * pheno

        # weights have 0 dimension if missing
        if(is.null(weights)) weights <- numeric(0)

        # null fit
        fit0 <- fit1_hk_addcovar(X0, # plug addcovar where genoprobs would be
                                 pheno,
                                 matrix(nrow=length(pheno), ncol=0),     # empty slot for addcovar
                                 weights, se=FALSE, tol)

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            fitA <- fit1_hk_addcovar(genoprobs, pheno, addcovar, weights, se=se, tol)
        }
        else {                  # intcovar
            fitA <- fit1_hk_intcovar(genoprobs, pheno, addcovar, intcovar,
                                     weights, se=se, tol)
        }

        # lod score
        n <- length(pheno)
        lod <- (n/2)*log10(fit0$rss/fitA$rss)

        # residual SDs using 1/n
        sigsq0 <- fit0$sigma^2/n*fit0$df
        sigsqA <- fitA$sigma^2/n*fitA$df

        # individual contributions to the lod score
        ind_lod <- 0.5*(fit0$resid^2/sigsq0 - fitA$resid^2/sigsqA + log(sigsq0) - log(sigsqA))/log(10)
        names(ind_lod) <- names(pheno)

        fitted <- pheno - fitA$resid
    }
    else { # binary phenotype
        # null fit
        if(is.null(weights)) { # no weights
            fit0 <- fit_binreg(X0, pheno, FALSE, maxit, bintol, tol, eta_max)
            weights <- numeric(0)
        }
        else {
            fit0 <- fit_binreg_weighted(X0, pheno, weights, FALSE, maxit, bintol, tol, eta_max)
        }

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            fitA <- fit1_binary_addcovar(genoprobs, pheno, addcovar, weights, se=se,
                                         maxit, bintol, tol, eta_max)
        }
        else {                  # intcovar
            fitA <- fit1_binary_intcovar(genoprobs, pheno, addcovar, intcovar,
                                         weights, se=se, maxit, bintol, tol, eta_max)
        }

        lod <- fitA$log10lik - fit0$log10lik

        # individual contributions to LOD
        if(length(weights)==0) weights <- rep(1, length(pheno))
        pA <- fitA$fitted_probs
        ind_lodA <- weights * (pheno*log10(pA) + (1-pheno)*log10(1-pA))
        p0 <- fit0$fitted_probs
        ind_lod0 <- weights * (pheno*log10(p0) + (1-pheno)*log10(1-p0))
        ind_lod <- ind_lodA - ind_lod0

        fitted <- stats::setNames(pA, names(pheno))
    }

    # names of coefficients
    coef_names <- scan1coef_names(genoprobs, addcovar, intcovar)

    # center the QTL effects at zero and add an intercept
    if(zerosum && is.null(contrasts)) {
        ng <- dim(genoprobs)[2]
        whval <- seq_len(ng)
        mu <- mean(fitA$coef[whval], na.rm=TRUE)
        fitA$coef <- c(fitA$coef, mu)
        fitA$coef[whval] <- fitA$coef[whval] - mu

        coef_names <- c(coef_names, "intercept")

        if(se) {
            fitA$SE <- c(fitA$SE, sqrt(mean(fitA$SE[whval]^2, na.rm=TRUE)))
        }
    }

    if(se) # results include standard errors
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names),
                    SE=stats::setNames(fitA$SE, coef_names),
                    fitted=fitted))
    else
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names),
                    fitted=fitted))
}
