# subset viterbi objects

#' Subsetting Viterbi results
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of \code{\link{viterbi}}
#'
#' @param x Imputed genotypes as output from \code{\link{viterbi}}.
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: logical values, or character
#' string IDs. Numbers are interpreted as character string IDs.
#' @param ... Ignored.
#'
#' @return The input viterbi object, with the selected
#' individuals and/or chromsomes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' \dontshow{grav2 <- grav2[1:8,c(1,2)]}
#' g <- viterbi(grav2)
#' # keep just individuals 1:5, chromosome 2
#' gsub <- g[1:5,2]
#' # keep just chromosome 2
#' gsub2 <- g[,2]
subset.viterbi <-
    function(x, ind=NULL, chr=NULL, ...)
{
    if(is.null(ind) && is.null(chr))
        stop("You must specify either ind or chr.")

    x_attr <- attributes(x)
    x_attrnam <- names(x_attr)
    x_class <- class(x)

    chrID <- names(x)
    n_chr <- length(chrID)
    if(!is.null(chr)) {
        if(is.logical(chr)) {
            if(length(chr) != n_chr) {
                stop("length(chr) [", length(chr), "] != no. chr in x [",
                     n_chr, "]")
            }
            chr <- chrID[chr] # convert to character strings
        }
        else {
            chrindex <- match(chr, chrID)
            if(any(is.na(chrindex))) {
                stop("Some chr not found: ",
                     paste(chr[is.na(chrindex)], collapse=", "))
            }
            chr <- chrID[chrindex] # character strings
        }
        if(length(chr) == 0)
            stop("Must retain at least one chromosome.")

        x <- unclass(x)[chr]
        to_sub <- "is_x_chr"
        for(obj in to_sub) {
            if(obj %in% x_attrnam)
                x_attr[[obj]] <- x_attr[[obj]][chr]
        }

    }

    indID <- rownames(x[[1]])

    n_ind <- length(indID)
    if(!is.null(ind)) {
        if(is.logical(ind)) {
            if(length(ind) != n_ind) {
                stop("length(ind) [", length(ind), "] != no. ind in x [",
                     n_ind, "]")
            }
            ind <- indID[ind] # convert to character strings
        }
        else if(is.numeric(ind)) {
            if(any(ind < 1 || ind > n_ind))
                stop("Numeric ind out of allowed range [1 - ", n_ind, "]")
            ind <- indID[ind] # convert to character strings
        }
        else {
            indindex <- match(ind, indID)
            if(any(is.na(indindex))) {
                stop("Some individuals not found: ",
                     paste(ind[is.na(indindex)], collapse=", "))
            }
            ind <- indID[indindex]
        }
        if(length(ind) == 0)
            stop("Must retain at least one individual.")

        for(i in names(x)) # loop over chromosomes
            x[[i]] <- unclass(x)[[i]][ind,,drop=FALSE]
    }

    for(obj in c("alleles", "is_x_chr", "crosstype"))
        attr(x, obj) <- x_attr[[obj]]
    class(x) <- x_class

    x
}

#' @export
#' @rdname subset.viterbi
`[.viterbi` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)
