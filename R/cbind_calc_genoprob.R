#' Join genotype probabilities for different chromosomes
#'
#' Join multiple genotype probability objects, as produced by
#' \code{\link{calc_genoprob}} for different individuals.
#'
#' @param ... Genotype probability objects as produced by
#' \code{\link{calc_genoprob}}. Must have the same set of individuals.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probsA <- calc_genoprob(grav2[1:5,1:2], step=1, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[1:5,3:4], step=1, error_prob=0.002)
#' probs <- cbind(probsA, probsB)
#'
#' @export
cbind.calc_genoprob <-
    function(...)
{
    args <- list(...)

    # to cbind: probs, map, grid, is_x_chr
    # to pass through (must match): is_female, cross_info, crosstype, alleles, alleleprobs, step, off_end, stepwidth

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    other_stuff <- c("is_female", "cross_info", "crosstype", "alleles", "alleleprobs",
                     "step", "off_end", "stepwidth", "error_prob", "map_function")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(args[[1]][[obj]], args[[i]][[obj]]))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
    }

    # paste stuff together
    main_stuff <- c("probs", "draws")
    for(i in 2:length(args)) {
        for(obj in c("probs", "draws")) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not present in all inputs")
            if(nrow(args[[1]][[obj]][[1]]) != nrow(args[[i]][[obj]][[1]]) ||
               !all(rownames(args[[1]][[obj]][[1]]) == rownames(args[[i]][[obj]][[1]])))
                stop("Input objects 1 and ", i, " have different individuals")

            result[[obj]] <- c(result[[obj]], args[[i]][[obj]])
        }
    }

    other_stuff <- c("map", "grid", "is_x_chr")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not present in all inputs")
            result[[obj]] <- c(result[[obj]], args[[i]][[obj]])
        }
    }

    result
}
