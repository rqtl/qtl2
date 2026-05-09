#' Join phased geno results for different individuals
#'
#' Join phased genotype objects, as produced by [guess_phase()],
#' for the same set of markers but for different individuals.
#'
#' @param ... Imputed genotype objects as produced by
#' [guess_phase()]. Must have the same set of markers.
#'
#' @return An object of class `"phasedgeno"`, like the input; see [guess_phase()].
#'
#' @seealso [cbind.phasedgeno()], [subset.phasedgeno()], [guess_phase()]
#'
#' @examples
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' m <- maxmarg(pr)
#' phA <- guess_phase(DOex[1:10,], m[1:10,])
#' phB <- guess_phase(DOex[11:20,], m[11:20,])
#' ph <- rbind(phA, phB)
#' }
#'
#' @export
rbind.phasedgeno <-
    rbind.calc_genoprob
