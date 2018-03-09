#' Join Viterbi results for different individuals
#'
#' Join multiple imputed genotype objects, as produced by
#' [viterbi()] for different individuals.
#'
#' @md
#'
#' @param ... Imputed genotype objects as produced by
#' [viterbi()]. Must have the same set of markers.
#'
#' @return A single viterbi object.
#'
#' @seealso [cbind.viterbi()], [viterbi()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' gA <- viterbi(grav2[1:5,], map, error_prob=0.002)
#' gB <- viterbi(grav2[6:12,], map, error_prob=0.002)
#' g <- rbind(gA, gB)
#'
#' @export
rbind.viterbi <-
    function(...)
{
    args <- list(...)

    # to rbind: the data itself
    # to pass through (must match): crosstype, is_x_chr, alleles

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    other_stuff <- c("crosstype", "is_x_chr", "alleles")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
    }

    for(i in 2:length(args)) {
        if(length(result) != length(args[[i]]) ||
           !all(names(result) == names(args[[i]])))
            stop("Input arguments have different chromosomes")

        for(s in seq(along=result)) {
            if(!is_same(ncol(result[[s]]), ncol(args[[i]][[s]])))
                stop("input objects have varying numbers of geno columns on chr ", names(result)[s])
            result[[s]] <- rbind(result[[s]], args[[i]][[s]])
        }
    }

    result
}
