# subset phasedgeno objects

#' Subsetting phased genotype objects
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of [guess_phase()]
#'
#' @param x Imputed genotypes as output from [guess_phase()], a list of 3d arrays.
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: logical values, or character
#' string IDs. Numbers are interpreted as character string IDs.
#' @param ... Ignored.
#'
#' @return An object of class `"phasedgeno"`, like the input, with the
#' selected individuals and/or chromosomes; see [guess_phase()].
#'
#' @seealso [guess_phase()], [cbind.phasedgeno()], [rbind.phasedgeno()]
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' m <- maxmarg(pr)
#' ph <- guess_phase(DOex, m)
#' ph <- ph[1:10,] # first ten individuals
#' ph <- ph[,2]  # chromosome 2
#' }
subset.phasedgeno <-
    subset.calc_genoprob

#' @export
#' @rdname subset.phasedgeno
`[.phasedgeno` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)
