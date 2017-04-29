#' Fit single-QTL model at a single position
#'
#' Fit a single-QTL model at a single putative QTL position and get detailed results
#' about estimated coefficients and individuals contributions to the LOD score.
#'
#' @param genoprobs A matrix of genotype probabilities, individuals x genotypes
#' @param pheno A numeric vector of phenotype values (just one phenotype, not a matrix of them)
#' @param kinship Optional kinship matrix.
#' @param addcovar An optional matrix of additive covariates.
#' @param nullcovar An optional matrix of additional additive
#' covariates that are used under the null hypothesis (of no QTL)
#' but not under the alternative (with a QTL). This is needed for
#' the X chromosome, where we might need sex as a additive
#' covariate under the null hypothesis, but we wouldn't want to
#' include it under the alternative as it would be collinear with
#' the QTL effects.
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have \code{names}
#' for individual identifiers. Ignored if \code{kinship} is provided.
#' @param contrasts An optional matrix of genotype contrasts, size
#' genotypes x genotypes. For an intercross, you might use
#' \code{cbind(c(1,1,1), c(-0.5, 0, 0.5), c(-0.5, 1, -0.5))} to get
#' mean, additive effect, and dominance effect. The default is the
#' identity matrix.
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If \code{model="binary"}, the phenotypes must have values in [0, 1].
#' @param se If TRUE, calculate the standard errors.
#' @param hsq (Optional) residual heritability; used only if
#' \code{kinship} provided.
#' @param reml If \code{kinship} provided: if \code{reml=TRUE}, use
#' REML; otherwise maximum likelihood.
#' @param tol Tolerance value for
#' linear regression by QR decomposition (in determining whether
#' columns are linearly dependent on others and should be omitted)
#' @param maxit Maximum number of iterations in logistic regression
#'     fit (when \code{model="binary"}).
#'
#' @return A list containing
#' \itemize{
#' \item \code{coef} - Vector of estimated coefficients.
#' \item \code{SE} - Vector of estimated standard errors (included if \code{se=TRUE}).
#' \item \code{lod} - The overall lod score.
#' \item \code{ind_lod} - Vector of individual contributions to the LOD score.
#' }
#'
#' @details For each of the inputs, the row names are used as
#' individual identifiers, to align individuals.
#'
#' If \code{kinship} is absent, Haley-Knott regression is performed.
#' If \code{kinship} is provided, a linear mixed model is used, with a
#' polygenic effect estimated under the null hypothesis of no (major)
#' QTL, and then taken as fixed as known in the genome scan.
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
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # leave-one-chromosome-out kinship matrix for chr 7
#' kinship7 <- calc_kinship(probs, "loco")[[7]]
#'
#' # scan chromosome 7 to find peak
#' out <- scan1(probs[,7], pheno, kinship7, addcovar=covar)
#'
#' # find peak position
#' max_pos <- rownames(max(out, map[7]))
#'
#' # fit QTL model just at that position
#' out_fit1 <- fit1(probs[[7]][,,max_pos], pheno, addcovar=covar)
#'
#' # fit QTL model just at that position, with polygenic effect
#' out_fit1_pg <- fit1(probs[[7]][,,max_pos], pheno, kinship7, addcovar=covar)
#'
#' @importFrom stats setNames
#' @export
fit1 <-
    function(genoprobs, pheno, kinship=NULL, addcovar=NULL, nullcovar=NULL,
             intcovar=NULL, weights=NULL,
             contrasts=NULL, model=c("normal", "binary"),
             se=TRUE, hsq=NULL, reml=TRUE, tol=1e-12, maxit=100)
{
    if(!is.null(kinship)) { # use LMM; see fit1_pg.R
        return(fit1_pg(genoprobs, pheno, kinship, addcovar, nullcovar,
                       intcovar, contrasts, se, hsq, reml, tol))
    }

    stopifnot(tol > 0)
    bintol <- sqrt(tol)

    # check that the objects have rownames
    check4names(pheno, addcovar, NULL, intcovar, nullcovar)

    # force things to be matrices
    if(!is.null(addcovar) && !is.matrix(addcovar))
        addcovar <- as.matrix(addcovar)
    if(!is.null(nullcovar) && !is.matrix(nullcovar))
        nullcovar <- as.matrix(nullcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
        intcovar <- as.matrix(intcovar)
    if(!is.null(contrasts) && !is.matrix(contrasts))
        contrasts <- as.matrix(contrasts)

    model <- match.arg(model)
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
    }
    else { # binary phenotype
        # null fit
        if(is.null(weights)) { # no weights
            fit0 <- fit_binreg(X0, pheno, FALSE, maxit, bintol, tol)
            weights <- numeric(0)
        }
        else {
            fit0 <- fit_binreg_weighted(X0, pheno, weights, FALSE, maxit, bintol, tol)
        }

        if(is.null(intcovar)) { # just addcovar
            if(is.null(addcovar)) addcovar <- matrix(nrow=length(ind2keep), ncol=0)
            fitA <- fit1_binary_addcovar(genoprobs, pheno, addcovar, weights, se=se,
                                         maxit, bintol, tol)
        }
        else {                  # intcovar
            fitA <- fit1_binary_intcovar(genoprobs, pheno, addcovar, intcovar,
                                         weights, se=se, maxit, bintol, tol)
        }

        lod <- fitA$log10lik - fit0$log10lik

        # individual contributions to LOD
        if(length(weights)==0) weights <- rep(1, length(pheno))
        pA <- fitA$fitted_probs
        ind_lodA <- weights * (pheno*log10(pA) + (1-pheno)*log10(1-pA))
        p0 <- fit0$fitted_probs
        ind_lod0 <- weights * (pheno*log10(p0) + (1-pheno)*log10(1-p0))
        ind_lod <- ind_lodA - ind_lod0
    }

    # names of coefficients
    coef_names <- scan1coef_names(genoprobs, addcovar, intcovar)

    if(se) # results include standard errors
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names),
                    SE=stats::setNames(fitA$SE, coef_names)))
    else
        return(list(lod=lod, ind_lod=ind_lod,
                    coef=stats::setNames(fitA$coef, coef_names)))
}
