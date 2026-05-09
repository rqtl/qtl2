#' Join phased genotype results for different chromosomes
#'
#' Join multiple phased genotype objects, as produced by [guess_phase()], for the
#' same set of individuals but different chromosomes.
#'
#' @param ... Imputed genotype objects as produced by
#' [guess_phase()]. Must have the same set of individuals.
#'
#' @return An object of class `"phasedgeno"`, like the input; see [guess_phase()].
#'
#' @seealso [rbind.phasedgeno()], [subset.phasedgeno()], [guess_phase()]
#'
#' @examples
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' m <- maxmarg(pr)
#' phA <- guess_phase(DOex[,2], m[,2])
#' phB <- guess_phase(DOex[,3], m[,3])
#' ph <- cbind(phA, phB)
#' }
#'
#' @export
cbind.phasedgeno <-
    cbind.calc_genoprob
