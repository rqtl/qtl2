#' Summarize scan1perm results
#'
#' Summarize permutation test results from [scan1perm()], as significance thresholds.
#'
#' @md
#'
#' @param object Output of [scan1perm()]
#' @param alpha Vector of significance levels
#'
#' @return
#' An object of class `summary.scan1perm`. If
#' [scan1perm()] was run with `perm_Xsp=FALSE`, this is
#' a single matrix of significance thresholds, with rows being
#' signicance levels and columns being the columns in the input. If
#' [scan1perm()] was run with `perm_Xsp=TRUE`, this is
#' a list of two matrices, with the significance thresholds for the
#' autosomes and X chromosome, respectively.
#'
#' The result has an attribute `"n_perm"` that has the numbers of
#' permutation replicates (either a matrix or a list of two matrices).
#'
#' @details
#' In the case of X-chromosome-specific permutations (when
#' [scan1perm()] was run with `perm_Xsp=TRUE`, we
#' follow the approach of Broman et al. (2006) to get separate
#' thresholds for the autosomes and X chromosome, using
#'
#' Let \eqn{L_A} and \eqn{L_X} be total the genetic lengths of the
#' autosomes and X chromosome, respectively, and let \eqn{L_T = L_A +
#' L_X} Then in place of \eqn{\alpha}{alpha}, we use \deqn{\alpha_A =
#' 1 - (1-\alpha)^{L_A/L_T}}{alpha_A = 1 - (1 - alpha)^(L_A/L_T)} as
#' the significance level for the autosomes and \deqn{\alpha_X = 1 -
#' (1-\alpha)^{L_X/L_T}}{alpha_x = 1 - (1 - alpha)^(LX/LT)} as the
#' significance level for the X chromosome.
#'
#' @references
#' Broman KW, Sen Ś, Owens SE, Manichaikul A, Southard-Smith EM,
#' Churchill GA (2006) The X chromosome in quantitative trait locus
#' mapping. Genetics 174:2151-2158
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c(10,18,"X")]}
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
#' # permutations with genome scan
#' \dontrun{
#' operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
#'                    n_perm=1000, perm_Xsp=TRUE,
#'                    chr_lengths=chr_lengths(iron$gmap))}
#' \dontshow{operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar, n_perm=3)}
#'
#' summary(operm, alpha=c(0.20, 0.05))
#'
#' @importFrom stats quantile
#' @export
summary_scan1perm <-
    function(object, alpha=0.05)
{

    if(is.matrix(object)) { # not X-chr-specific
        result <- apply(object, 2, quantile, 1-alpha, na.rm=TRUE)
        if(length(alpha)==1) {
            result <- matrix(result, nrow=1)
            colnames(result) <- colnames(object)
        }

        rownames(result) <- alpha

        n_perm <- matrix(colSums(!is.na(object)), nrow=1)
        colnames(n_perm) <- colnames(object)

        class(result) <- c("summary.scan1perm", "matrix")
    } else { # X-chr specific
        if(!is.list(object) || length(object) != 2 || any(sort(names(object)) != c("A", "X")))
            stop('input should be a matrix or a list of two matrices named "A" and "X"')

        chr_lengths <- attr(object, "chr_lengths")
        if(is.null(chr_lengths) || length(chr_lengths) != 2
           || any(names(chr_lengths) != c("A", "X"))) {
            warning("Ill-formed chr_lengths attribute.\n",
                    "Estimating LA/LX from the numbers of permutations, as ",
                    round(nrow(object$A)/nrow(object$X), 2))
            LA <- nrow(object$A)
            LX <- nrow(object$X)
        } else {
            LA <- chr_lengths["A"]
            LX <- chr_lengths["X"]
        }
        Lt <- LA+LX

        alphaA <- 1 - (1-alpha)^(LA/Lt)
        alphaX <- 1 - (1-alpha)^(LX/Lt)

        result <- list(A=apply(object$A, 2, quantile, 1-alphaA, na.rm=TRUE),
                       X=apply(object$X, 2, quantile, 1-alphaX, na.rm=TRUE))
        for(i in seq_along(result)) {
            if(!is.matrix(result[[i]])) {
                result[[i]] <- matrix(result[[i]], nrow=1)
                colnames(result[[i]]) <- colnames(object[[i]])
            }
            rownames(result[[i]]) <- alpha
        }
        n_perm <- rbind(A=colSums(!is.na(object$A)),
                        X=colSums(!is.na(object$X)))
        colnames(n_perm) <- colnames(object$A)

        class(result) <- c("summary.scan1perm", "list")
    }

    attr(result, "n_perm") <- n_perm
    result
}

#' @rdname summary_scan1perm
#' @param ... Ignored
#' @export
summary.scan1perm <-
    function(object, alpha=0.05, ...)
{
    summary_scan1perm(object, alpha=alpha)
}

#' Print summary of scan1perm permutations
#'
#' Print summary of scan1perm permutations
#'
#' @md
#'
#' @param x Object of class `"summary.scan1perm"`, as produced by [summary_scan1perm()].
#' @param digits Number of digits in printing significance thresholds; passed to [base::print()].
#' @param ... Ignored.
#'
#' @return Invisibly returns the input, `x`.
#'
#' @details This is to go with [summary_scan1perm()], so
#' that the summary output is printed in a nice format. Generally
#' not called directly, but it can be in order to control the
#' number of digits that appear.
#'
#' @examples
#' \donttest{
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
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
#' # permutations with genome scan
#' operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
#'                    n_perm=100, perm_Xsp=TRUE,
#'                    chr_lengths=chr_lengths(iron$gmap))
#'
#' summary(operm, alpha=c(0.20, 0.05))
#'
#' print( summary(operm, alpha=c(0.20, 0.05)), digits=8 )
#' }
#'
#' @export
print.summary.scan1perm <-
    function(x, digits=3, ...)
{
    n_perm <- attr(x, "n_perm")

    if(is.list(x)) { # X-chr-specific
        if(length(unique(n_perm[1,])) == 1 &&
           length(unique(n_perm[2,])) == 1) constant_perms <- TRUE
        else constant_perms <- FALSE

        cat("Autosome LOD thresholds ")
        if(constant_perms) cat("(", n_perm[1,1], " permutations)", sep="")
        cat("\n")
        print(x$A, digits=digits)

        cat("\nX chromosome LOD thresholds ")
        if(constant_perms) cat("(", n_perm[2,1], " permutations)", sep="")
        cat("\n")
        print(x$X, digits=digits)

        if(!constant_perms) {
            cat("\n")
            cat("Number of permutations:\n")
            n_spaces <- max(vapply(x, function(a) nchar(rownames(a)), 1))
            rownames(n_perm) <- paste0(paste(rep(" ", n_spaces-1), collapse=""), rownames(n_perm))
            print(n_perm)
        }
    }
    else { # not X-chr-specific
        if(length(unique(as.numeric(n_perm))) == 1) constant_perms <- TRUE
        else constant_perms <- FALSE

        cat("LOD thresholds ")
        if(constant_perms) cat("(", n_perm[1], " permutations)", sep="")
        cat("\n")

        print(x[seq_len(nrow(x)),,drop=FALSE], digits=digits)

        if(!constant_perms) {
            cat("\n")
            cat("Number of permutations:\n")
            n_spaces <- max(nchar(rownames(x)))
            rownames(n_perm) <- paste(rep(" ", n_spaces), collapse="")
            print(n_perm)
        }
    }

    invisible(x)
}
