#' Join viterbi results for different chromosomes
#'
#' Join multiple viterbi objects, as produced by [viterbi()], for the
#' same set of individuals but different chromosomes.
#'
#' @param ... Imputed genotype objects as produced by
#' [viterbi()]. Must have the same set of individuals.
#'
#' @return An object of class `"viterbi"`, like the input; see [viterbi()].
#'
#' @seealso [rbind.viterbi()], [viterbi()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' gA <- viterbi(grav2[1:5,1:2], map, error_prob=0.002)
#' gB <- viterbi(grav2[1:5,3:4], map, error_prob=0.002)
#' g <- cbind(gA, gB)
#'
#' @export
cbind.viterbi <-
    cbind.calc_genoprob
