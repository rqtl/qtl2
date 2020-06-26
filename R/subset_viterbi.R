# subset viterbi objects

#' Subsetting Viterbi results
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of [viterbi()]
#'
#' @param x Imputed genotypes as output from [viterbi()].
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: logical values, or character
#' string IDs. Numbers are interpreted as character string IDs.
#' @param ... Ignored.
#'
#' @return An object of class `"viterbi"`, like the input, with the
#' selected individuals and/or chromosomes; see [viterbi()].
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
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

    if(!is.null(chr)) {
        chr <- subset_chr(chr, names(x))
        if(length(chr) == 0)
            stop("Must retain at least one chromosome.")

        x <- unclass(x)[chr]
        to_sub <- "is_x_chr"
        for(obj in to_sub) {
            if(obj %in% x_attrnam)
                x_attr[[obj]] <- x_attr[[obj]][chr]
        }

    }

    if(!is.null(ind)) {
        ind <- subset_ind(ind, rownames(x[[1]]))
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
