#' Combine data from scan1perm objects
#'
#' Row-bind multiple scan1perm objects with the same set of columns
#'
#' @param ... A set of permutation results from
#' \code{\link{scan1perm}} (objects of class \code{scan1perm}).
#' They must have the same set of columns. If any include
#' autosome/X chromosome-specific permutations, they must all be
#' such.
#'
#' @return The combined row-binded input, as a \code{scan1perm} object.
#'
#' @details The aim of this function is to concatenate the results
#' from multiple runs of a permutation test with
#' \code{\link{scan1perm}}, to assist in the case that such
#' permutations are done on multiple processors in parallel.
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
#' # permutations with genome scan
#' \dontrun{
#' operm1 <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
#'                 n_perm=500, perm_Xsp=TRUE,
#'                 chr_lengths=chr_lengths(iron$gmap))
#' operm2 <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
#'                 n_perm=500, perm_Xsp=TRUE,
#'                 chr_lengths=chr_lengths(iron$gmap))
#' }
#'
#' \dontshow{
#' operm1 <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar, n_perm=3)
#' operm2 <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar, n_perm=3)
#' }
#'
#' operm <- rbind(operm1, operm2)
#'
#' @seealso \code{\link{cbind.scan1perm}}, \code{\link{scan1perm}}, \code{\link{scan1}}
#'
#' @export
rbind.scan1perm <-
    function(...)
{
    dots <- list(...)

    if(length(dots)==1) return(dots[[1]])

    if(is.matrix(dots[[1]])) {
        if(!all(vapply(dots, is.matrix, TRUE)))
            stop("Inputs cannot be a mixture of X-chr-specific permutations and not")
        n_col <- vapply(dots, ncol, 1)
        if(length(unique(n_col)) != 1)
            stop("Inputs must all have the same number of columns.")
        for(i in seq_along(dots)[-1]) {
            if(!all(colnames(dots[[1]]) == colnames(dots[[i]])))
                warning("Inputs have different column names; using those from the first input.")
        }

        result <- do.call("rbind", lapply(dots, unclass))
        colnames(result) <- colnames(dots[[1]])
        class(result) <- class(dots[[1]])
        return(result)
    }
    else { # perm_Xsp=TRUE
        # check they're all lists of A and X
        if(!all(vapply(dots, is.list, TRUE)))
            stop("Inputs cannot be a mixture of X-chr-specific permutations and not")
        if(!all(vapply(dots, length, 1)==2) || !all(vapply(dots, function(a) all(sort(names(a)) == c("A", "X")), TRUE)))
            stop("Inputs should all be a list with two components")

        # check that the columns are the same
        n_colA <- vapply(dots, function(a) ncol(a$A), 1)
        n_colX <- vapply(dots, function(a) ncol(a$X), 1)
        if(length(unique(c(n_colA, n_colX))) != 1)
            stop("Inputs must all have the same number of columns.")
        for(i in seq_along(dots)[-1]) {
            if(!all(colnames(dots[[1]]$A) == colnames(dots[[i]]$A)) ||
               !all(colnames(dots[[1]]$X) == colnames(dots[[i]]$X)))
                warning("Inputs have different column names; using those from the first input.")
        }

        # check the chr_lengths attributes
        chr_lengths <- lapply(dots, function(a) attr(a, "chr_lengths"))
        if(any(vapply(chr_lengths, is.null, TRUE)))
            stop("Missing chr_lengths attribute from some inputs")
        if(any(vapply(chr_lengths, length, 1) != 2) ||
           any(vapply(chr_lengths, function(a) any(sort(names(a)) != c("A", "X")), TRUE)))
            stop('chr_lengths attributes should be length 2 vectors with names "A" and "X"')
        if(max(vapply(chr_lengths, function(a) max(abs(chr_lengths[[1]] - a)), 1.0)) > 1e-6) {
            stop("chr_lengths attributes differ.")
        }

        # check the is_x_chr attributes of the chr_lengths attributes
        is_x_chr <- lapply(chr_lengths, function(a) attr(a, "is_x_chr"))
        if(!any(vapply(is_x_chr, is.null, TRUE))) { # if some missing; just skip this
            if(any(vapply(is_x_chr, length, 1) != 2) ||
               any(vapply(is_x_chr, function(a) any(sort(names(a)) != c("A", "X")), TRUE)))
                stop('is_x_chr attributes should be length 2 vectors with names "A" and "X"')
            for(i in seq_along(is_x_chr)[-1]) {
                if(is_x_chr[[1]]["A"] != is_x_chr[[i]]["A"] ||
                   is_x_chr[[1]]["X"] != is_x_chr[[i]]["X"])
                    stop("is_x_chr attributes differ")
            }
            is_x_chr <- is_x_chr[[1]]
        }
        else is_x_chr <- NULL # if some missing, just skip it

        # rbind
        A <- do.call("rbind", lapply(dots, function(a) a$A))
        X <- do.call("rbind", lapply(dots, function(a) a$X))

        # add attributes
        colnames(A) <- colnames(dots[[1]]$A)
        colnames(X) <- colnames(dots[[1]]$X)
        result <- list(A=A, X=X)
        attr(chr_lengths[[1]], "is_x_chr") <- is_x_chr
        attr(result, "chr_lengths") <- chr_lengths[[1]]
        class(result) <- class(dots[[1]])

        return(result)
    }

}


#' @export
#' @rdname rbind.scan1perm
c.scan1perm <- rbind.scan1perm


#' Combine columns from multiple scan1 permutation results
#'
#' Column-bind multiple scan1perm objects with the same numbers of rows.
#'
#' @param ... A set of permutation results from
#' \code{\link{scan1perm}} (objects of class \code{scan1perm}). If
#' different numbers of permutation replicates were used, those
#' columns with fewer replicates are padded with missing values
#' \code{NA}. However, if any include autosome/X
#' chromosome-specific permutations, they must all be such.
#'
#' @return The combined column-binded input, as a \code{scan1perm} object.
#'
#' @details The aim of this function is to concatenate the results
#' from multiple runs of a permutation test with
#' \code{\link{scan1perm}}, generally with different phenotypes
#' and/or methods, to be used in parallel with
#' \code{\link{rbind.scan1perm}}.
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
#' # permutations with genome scan
#' \dontrun{
#' operm1 <- scan1perm(probs, pheno[,1,drop=FALSE], addcovar=covar, Xcovar=Xcovar,
#'                 n_perm=1000, perm_Xsp=TRUE,
#'                 chr_lengths=chr_lengths(iron$gmap))
#' operm2 <- scan1perm(probs, pheno[,2,drop=FALSE], addcovar=covar, Xcovar=Xcovar,
#'                 n_perm=1000, perm_Xsp=TRUE,
#'                 chr_lengths=chr_lengths(iron$gmap))
#' }
#'
#' \dontshow{
#' operm1 <- scan1perm(probs, pheno[,1,drop=FALSE], addcovar=covar, Xcovar=Xcovar, n_perm=3)
#' operm2 <- scan1perm(probs, pheno[,2,drop=FALSE], addcovar=covar, Xcovar=Xcovar, n_perm=3)
#' }
#'
#' operm <- cbind(operm1, operm2)
#'
#' @seealso \code{\link{rbind.scan1perm}}, \code{\link{scan1perm}}, \code{\link{scan1}}
#'
#' @export
cbind.scan1perm <-
    function(...)
{
    dots <- list(...)

    if(length(dots)==1) return(dots[[1]])

    if(is.matrix(dots[[1]])) {
        if(!all(vapply(dots, is.matrix, TRUE)))
            stop("Inputs cannot be a mixture of X-chr-specific permutations and not")

        result <- do.call("cbind_expand", lapply(dots, unclass))
        class(result) <- class(dots[[1]])
        return(result)
    }
    else { # perm_Xsp=TRUE
        # check they're all lists of A and X
        if(!all(vapply(dots, is.list, TRUE)))
            stop("Inputs cannot be a mixture of X-chr-specific permutations and not")
        if(!all(vapply(dots, length, 1)==2) || !all(vapply(dots, function(a) all(sort(names(a)) == c("A", "X")), TRUE)))
            stop("Inputs should all be a list with two components")

        # check the chr_lengths attributes
        chr_lengths <- lapply(dots, function(a) attr(a, "chr_lengths"))
        if(any(vapply(chr_lengths, is.null, TRUE)))
            stop("Missing chr_lengths attribute from some inputs")
        if(any(vapply(chr_lengths, length, 1) != 2) ||
           any(vapply(chr_lengths, function(a) any(sort(names(a)) != c("A", "X")), TRUE)))
            stop('chr_lengths attributes should be length 2 vectors with names "A" and "X"')
        if(max(vapply(chr_lengths, function(a) max(abs(chr_lengths[[1]] - a)), 1.0)) > 1e-6) {
            stop("chr_lengths attributes differ.")
        }

        # check the is_x_chr attributes of the chr_lengths attributes
        is_x_chr <- lapply(chr_lengths, function(a) attr(a, "is_x_chr"))
        if(!any(vapply(is_x_chr, is.null, TRUE))) { # if some missing; just skip this
            if(any(vapply(is_x_chr, length, 1) != 2) ||
               any(vapply(is_x_chr, function(a) any(sort(names(a)) != c("A", "X")), TRUE)))
                stop('is_x_chr attributes should be length 2 vectors with names "A" and "X"')
            for(i in seq_along(is_x_chr)[-1]) {
                if(is_x_chr[[1]]["A"] != is_x_chr[[i]]["A"] ||
                   is_x_chr[[1]]["X"] != is_x_chr[[i]]["X"])
                    stop("is_x_chr attributes differ")
            }
            is_x_chr <- is_x_chr[[1]]
        }
        else is_x_chr <- NULL # if some missing, just skip it

        # rbind
        A <- do.call("cbind_expand", lapply(dots, function(a) a$A))
        X <- do.call("cbind_expand", lapply(dots, function(a) a$X))

        # add attributes
        result <- list(A=A, X=X)
        attr(chr_lengths[[1]], "is_x_chr") <- is_x_chr
        attr(result, "chr_lengths") <- chr_lengths[[1]]
        class(result) <- class(dots[[1]])

        return(result)
    }

}
