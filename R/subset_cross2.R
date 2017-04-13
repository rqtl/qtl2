# subset cross2 objects

# subset.cross2
#' Subsetting data for a QTL experiment
#'
#' Pull out a specified set of individuals and/or chromosomes from a
#' \code{cross2} object.
#'
#' @param x An object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs.
#' @param chr A vector of chromosomes: numeric indices, logical
#' values, or character string IDs
#' @param ... Ignored.
#'
#' @details
#' When subsetting by individual, if \code{ind} is numeric, they're
#' assumed to be numeric indices; if character strings, they're
#' assumed to be individual IDs. \code{ind} can be numeric or logical
#' only if the genotype, phenotype, and covariate data all have the
#' same individuals in the same order.
#'
#' When subsetting by chromosome, \code{chr} is \emph{always}
#' converted to character strings and treated as chromosome IDs. So if
#' there are three chromosomes with IDs \code{"18"}, \code{"19"}, and
#' \code{"X"}, \code{mycross[,18]} will give the first of the
#' chromosomes (labeled \code{"18"}) and \code{mycross[,3]} will give
#' an error.
#'
#' When using character string IDs for \code{ind} or \code{chr}, you
#' can use "negative" subscripts to indicate exclusions, for example
#' \code{mycross[,c("-18", "-X")]} or \code{mycross["-Mouse2501",]}.
#' But you can't mix "positive" and "negative" subscripts, and if any
#' of the individuals has an ID that begins with \code{"-"}, you can't
#' use negative subscripts like this.
#'
#' @return The input \code{cross2} object, with the selected
#' individuals and/or chromsomes.
#'
#' @section Warning:
#' The order of the two arguments is reversed relative to the related
#' function in \href{http://rqtl.org}{R/qtl}.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' # keep individuals 1-20 and chromosomes 3 and 4
#' grav2sub <- grav2[1:20, c(3,4)]
#' # keep just chromosome 1
#' grav2_c1 <- grav2[,1]
subset.cross2 <-
    function(x, ind=NULL, chr=NULL, ...)
{
    if(is.null(ind) && is.null(chr))
        stop("You must specify either ind or chr.")

    slice_by_chr <- c("geno", "founder_geno", "gmap", "pmap", "is_x_chr")

    slice_by_ind <- c("geno", "is_female", "cross_info", "pheno", "covar")

    if(!is.null(chr)) {
        # first clean up chromosome argument
        chr <- subset_chr(chr, names(x$geno))

        for(obj in slice_by_chr) {
            if(obj %in% names(x)) {
                x[[obj]] <- x[[obj]][chr]
            }
        }
    }

    if(!is.null(ind)) {
        # check that geno, covar, pheno have individuals in same order
        allow_logical <- TRUE
        ind_g <- rownames(x$geno[[1]])
        if(!is.null(x$covar) && (nrow(x$covar) != length(ind_g) ||
           !all(rownames(x$covar) == ind_g)))
            allow_logical <- FALSE
        if(!is.null(x$pheno) && (nrow(x$pheno) != length(ind_g) ||
           !all(rownames(x$pheno) == ind_g)))
            allow_logical <- FALSE

        # first clean up ind argument
        if(is.logical(ind)) {
            if(!allow_logical) stop("ind can't be logical if different individuals in geno, pheno, covar")

            if(length(ind) != length(ind_g))
                stop("ind is logical but length [", length(ind), "] != no. ind in cross [",
                     length(ind_g), "]")
            ind <- ind_g[ind]
        }
        if(is.numeric(ind)) { # treat as numeric indexes
            if(!allow_logical) stop("ind can't be numeric if different individuals in geno, pheno, covar")
            if(any(ind < 0)) { # deal with negatives
                if(!all(ind < 0)) stop("Can't mix negative and positive ind subscripts")
                ind <- (seq_along(ind_g))[ind]
            }

            if(any(ind < 1 | ind > length(ind_g))) {
                ind <- ind[ind >= 1 & ind <= length(ind_g)]
                if(length(ind)==0) stop("All ind out of range")
                warning("some ind out of range [1,", length(ind_g), "]")
            }

            ind <- ind_g[ind]
        }

        # now treat ind as character strings
        ind <- as.character(ind)
        allind <- ind_ids(x)

        # look for negatives; turn to positives
        if(!any(grepl("^\\-", allind)) && any(grepl("^\\-", ind))) { # if "-" used in actual IDs, don't allow negative subscripts
            if(!all(grepl("^\\-", ind)))
                stop("Can't mix negative and positive ind subscripts")
            ind <- sub("^\\-", "", ind)
            if(!all(ind %in% allind)) {
                if(!any(ind %in% allind)) stop("None of the individuals in cross")
                warning("Some individuals not in cross: ", paste(ind[!(ind %in% allind)], collapse=", "))
            }
            ind <- allind[!(allind %in% ind)]
        }

        if(!all(ind %in% allind)) {
            if(!any(ind %in% allind)) stop("None of the individuals in cross")
            warning("Some ind not in cross: ", paste(ind[!(ind %in% allind)], collapse=", "))
            ind <- ind[ind %in% ind]
        }

        # Finally, the actual subsetting
        for(obj in slice_by_ind) {
            if(obj %in% names(x)) {
                # is it a matrix like cross_info?
                if(is.matrix(x[[obj]]) || is.data.frame(x[[obj]])) {
                    x[[obj]] <- x[[obj]][ind[ind %in% rownames(x[[obj]])],,drop=FALSE]
                }
                # is it a list like $geno?
                else if(!is.matrix(x[[obj]]) && is.list(x[[obj]])) {
                    x[[obj]] <- lapply(x[[obj]], function(a, b) {
                        if(is.matrix(a))
                            a <- a[b[b %in% rownames(a)],,drop=FALSE]
                        else
                            a <- a[b[b %in% names(a)]]
                        a }, ind)
                }
                else # it's a vector
                    x[[obj]] <- x[[obj]][ind[ind %in% names(x[[obj]])]]
            }
        }
    }

    x
}

#' @export
#' @rdname subset.cross2
`[.cross2` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)
