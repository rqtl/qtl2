#' Join viterbi results for different chromosomes
#'
#' Join multiple viterbi objects, as produced by
#' \code{\link{viterbi}} for different individuals.
#'
#' @param ... Imputed genotype objects as produced by
#' \code{\link{viterbi}}. Must have the same set of individuals.
#'
#' @return A single imputed genotype object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' gA <- viterbi(grav2[1:5,1:2], step=1, error_prob=0.002)
#' gB <- viterbi(grav2[1:5,3:4], step=1, error_prob=0.002)
#' g <- cbind(gA, gB)
#'
#' @export
cbind.viterbi <-
    function(...)
{
    args <- list(...)

    # to cbind: the data itself, is_x_chr
    # to pass through (must match): crosstype, alleles

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # paste stuff together
    for(i in 2:length(args)) {
        if(nrow(args[[1]][[1]]) != nrow(args[[i]][[1]]) ||
           !all(rownames(args[[1]][[1]]) == rownames(args[[i]][[1]])))
            stop("Input objects 1 and ", i, " have different individuals")

        result <- c(result, args[[i]])
    }

    other_stuff <- "is_x_chr"
    for(obj in other_stuff)
        attr(result, obj) <- attr(args[[1]], obj)
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(is.null(attr(args[[1]], obj)) && is.null(attr(args[[i]], obj))) next # not present
            if(is.null(attr(args[[1]], obj)) || is.null(attr(args[[i]], obj)))
                stop(obj, " not present in all inputs")
            attr(result, obj) <- c(attr(result, obj), attr(args[[i]], obj))
        }
    }

    # check that things match
    other_stuff <- c("crosstype", "alleles")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
    }
    for(obj in other_stuff)
        attr(result, obj) <- attr(args[[1]], obj)

    class(result) <- class(args[[1]])

    result
}
