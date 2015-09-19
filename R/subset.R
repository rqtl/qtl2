# subset cross2, calc_genoprob, and sim_geno objects

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
#' function in \href{http://www.rqtl.org}{R/qtl}.
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
    function(x, ind, chr, ...)
{
    if(missing(ind) && missing(chr))
        stop("You must specify either ind or chr.")

    slice_by_chr <- c("geno", "founder_geno", "gmap", "pmap", "is_x_chr")

    slice_by_ind <- c("geno", "is_female", "cross_info")
    slice_by_ind_linemap <- c("pheno", "covar")

    if(!missing(chr)) {
        for(obj in slice_by_chr) {
            if(obj %in% names(x)) {
                x[[obj]] <- x[[obj]][chr]
            }
        }
    }

    if(!missing(ind)) {
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
    function(x, ind, chr)
    subset(x, ind, chr)

# subset.calc_genoprob
#' Subsetting genotype probabilities
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of \code{\link{calc_genoprob}}.
#'
#' @param x A list of arrays of genotype probabilities, as output from \code{\link{calc_genoprob}}.
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: numeric indices, logical
#' values, or character string IDs
#' @param ... Ignored.
#'
#' @return The input genotype probabilities, with the selected
#' individuals and/or chromsomes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' \dontshow{grav2 <- grav2[1:8,c(1,2)]}
#' pr <- calc_genoprob(grav2)
#' # keep just individuals 1:5, chromosome 2
#' prsub <- pr[1:5,2]
#' # keep just chromosome 2
#' prsub2 <- pr[,2]
subset.calc_genoprob <-
    function(x, ind, chr, ...)
{
    if(missing(ind) && missing(chr))
        stop("You must specify either ind or chr.")

    at <- attributes(x)
    atnam <- names(at)

    if(!missing(chr)) {
        x <- unclass(x)[chr]

        to_sub <- c("map", "is_x_chr")
        for(objname in to_sub) {
            if(objname %in% atnam)
                at[[objname]] <- at[[objname]][chr]
        }
    }

    if(!missing(ind)) {
        x <- lapply(x, "[", ind, , , drop=FALSE)

        if("cross_info" %in% atnam)
            at$cross_info <- at$cross_info[ind, , drop=FALSE]
    }

    # restore attributes
    for(a in atnam[atnam != "names"])
        attr(x, a) <- at[[a]]

    x
}

#' @export
#' @rdname subset.calc_genoprob
`[.calc_genoprob` <-
    function(x, ind, chr)
    subset(x, ind, chr)

# subset.sim_geno
#' Subsetting imputed genotypes
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of \code{\link{sim_geno}}.
#'
#' @param x A list of arrays of imputed genotypes, as output from \code{\link{sim_geno}}
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: numeric indices, logical
#' values, or character string IDs
#' @param ... Ignored.
#'
#' @return The input imputed genotypes, with the selected
#' individuals and/or chromsomes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' \dontshow{grav2 <- grav2[1:8,c(1,2)]}
#' dr <- sim_geno(grav2, n_draws=4)
#' # keep just individuals 1:5, chromosome 2
#' drsub <- dr[1:5,2]
#' # keep just chromosome 2
#' drsub2 <- dr[,2]
subset.sim_geno <- subset.calc_genoprob

#' @export
#' @rdname subset.sim_geno
`[.sim_geno` <- `[.calc_genoprob`
