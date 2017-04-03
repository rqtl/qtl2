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
#' values, or character string IDs
#' @param chr A vector of chromosomes: numeric indices, logical
#' values, or character string IDs
#' @param ... Ignored.
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

    slice_by_ind <- c("geno", "is_female", "cross_info")
    slice_by_ind_linemap <- c("pheno", "covar")

    if(!is.null(chr)) {
        # first clean up chromosome argument
        allchr <- names(x$geno)
        if(is.logical(chr)) {
            if(length(chr) != length(allchr))
                stop("chr is logical but length [", length(chr), "] != n_chr is x [",
                     length(allchr), "]")
            chr <- allchr[chr]
        } else {
            chr <- as.character(chr)
            if(!all(chr %in% allchr)) {
                if(!any(chr %in% allchr)) stop("None of the chromosomes in x")
                warning("Some chr not in x: ", paste(chr[!(chr %in% allchr)], collapse=", "))
                chr <- chr[chr %in% allchr]
            }
        }

        for(obj in slice_by_chr) {
            if(obj %in% names(x)) {
                x[[obj]] <- x[[obj]][chr]
            }
        }
    }

    if(!is.null(ind)) {
        # first clean up ind argument
        allind <- rownames(x$geno[[1]])
        if(is.logical(ind)) {
            if(length(ind) != allind)
                stop("ind is logical but length [", length(ind), "] != n_ind is x [",
                     length(allind), "]")
            ind <- allind[ind]
        }
        if(is.numeric(ind)) { # treat as numeric indexes
            if(any(ind < 1 | ind > length(allind))) {
                ind <- ind[ind >= 1 & ind <= length(allind)]
                if(length(ind)==0)
                    stop("All ind out of range")
                warning("some ind out of range [1,", length(allind), "]")
            }
        }
        else { # character
            ind <- as.character(ind)
            if(!all(ind %in% allind)) {
                if(!any(ind %in% allind)) stop("None of the individuals in x")
                warning("Some ind not in x: ", paste(ind[!(ind %in% allind)], collapse=", "))
                ind <- ind[ind %in% allind]
            }
        }

        for(obj in slice_by_ind) {
            if(obj %in% names(x)) {
                # is it a list like $geno?
                if(!is.matrix(x[[obj]]) && is.list(x[[obj]])) {
                    x[[obj]] <- lapply(x[[obj]], function(a, b) {
                        if(is.matrix(a))
                            a <- a[b,,drop=FALSE]
                        else
                            a <- a[b]
                        a }, ind)
                }
                # is it a matrix like cross_info?
                else if(is.matrix(x[[obj]]))
                    x[[obj]] <- x[[obj]][ind,,drop=FALSE]
                else # it's a vector
                    x[[obj]] <- x[[obj]][ind]
            }
        }

        if("linemap" %in% names(x)) {
            # slice the lines then the individuals by lines
            ind2keep <- names(x$linemap)[x$linemap %in% rownames(x$geno)]
            x$linemap <- x$linemap[ind2keep]
            for(obj in slice_by_ind_linemap) {
                x[[obj]] <- x[[obj]][ind2keep,,drop=FALSE]
            }
        }
        else {
            for(obj in slice_by_ind_linemap) {
                x[[obj]] <- x[[obj]][ind,,drop=FALSE]
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
