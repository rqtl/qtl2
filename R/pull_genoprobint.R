#' Pull genotype probabilities for an interval
#'
#' Pull out the genotype probabilities for a given genomic interval
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param map The marker map for the genotype probabilities
#' @param chr Chromosome ID (single character sting)
#' @param interval Interval (pair of numbers)
#'
#' @return A list containing a single 3d array of genotype probabilities, like the input `genoprobs`
#' but for the designated interval.
#'
#' @seealso [find_marker()], [pull_genoprobpos()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,"8"]}
#' gmap <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, gmap, error_prob=0.002)
#'
#' pr_sub <- pull_genoprobint(pr, gmap, "8", c(25, 35))
#'
#' @export
pull_genoprobint <-
    function(genoprobs, map, chr, interval)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")

    # get vector of marker names that will be kept
    markers <- find_marker(map, chr, interval=interval)
    if(length(markers)==0) stop("No markers/pseudomarkers in the interval")

    # reduce to the one chromosome
    genoprobs <- genoprobs[,chr]

    # reduce to the markers in the interval
    genoprobs[[1]] <- genoprobs[[1]][,,markers,drop=FALSE]

    genoprobs
}
