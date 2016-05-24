#' Join Viterbi results for different individuals
#'
#' Join multiple imputed genotype objects, as produced by
#' \code{\link{viterbi}} for different individuals.
#'
#' @param ... Imputed genotype objects as produced by
#' \code{\link{viterbi}}. Must have the same set of markers.
#'
#' @return A single viterbi object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' gA <- viterbi(grav2[1:5,], step=1, error_prob=0.002)
#' gB <- viterbi(grav2[6:12,], step=1, error_prob=0.002)
#' g <- rbind(gA, gB)
#'
#' @export
rbind.viterbi <-
    function(...)
{
    args <- list(...)

    # to rbind: geno, sex, cross_info
    # to pass through (must match): map, grid, crosstype, is_x_chr, alleles, alleleprobs, step, off_end, stepwidth

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    nested_stuff <- c("map", "grid")
    other_stuff <- c("crosstype", "is_x_chr", "alleles", "alleleprobs",
                     "step", "off_end", "stepwidth", "error_prob", "map_function")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(args[[1]][[obj]], args[[i]][[obj]]))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
        for(obj in nested_stuff) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not prsent in all inputs")
            if(!is_same(names(args[[1]][[obj]]), names(args[[i]][[obj]])))
                stop("Input objects 1 and ", i, " differ in their ", obj)
            for(chr in seq(along=args[[1]][[obj]])) {
                if(!is_same(args[[1]][[obj]][[chr]], args[[i]][[obj]][[chr]]))
                    stop("Input objects 1 and ", i, " differ in their ", obj,
                         " on chromosome ", chr)
            }
        }
    }

    for(i in 2:length(args)) {
        if(!("sex" %in% names(result)) && !("sex" %in% names(args[[i]]))) next
        if(!("sex" %in% names(result) && "sex" %in% names(args[[i]])))
            stop("sex present in only some input objects")

        if(!("cross_info" %in% names(result)) && !("cross_info" %in% names(args[[i]]))) next
        if(!("cross_info" %in% names(result) && "cross_info" %in% names(args[[i]])))
            stop("cross_info present in only some input objects")
        if(!is_same(ncol(result$cross_info), ncol(args[[i]]$cross_info)))
            stop("input objects have varying numbers of cross_info columns")
        result$cross_info <- rbind(result$cross_info, args[[i]]$cross_info)

        if(!("geno" %in% names(result) && "geno" %in% names(args[[i]])))
            stop("geno present in only some input objects")
        for(s in seq(along=result$geno)) {
            if(!is_same(ncol(result$geno[[s]]), ncol(args[[i]]$geno[[s]])))
                stop("input objects have varying numbers of geno columns on chr ", names(result$geno)[s])
            result$geno[[s]] <- rbind(result$geno[[s]], args[[i]]$geno[[s]])
        }
    }

    result
}
