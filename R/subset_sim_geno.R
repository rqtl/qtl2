# subset sim_geno objects

#' Subsetting imputed genotypes
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of [sim_geno()].
#'
#' @md
#'
#' @param x Imputed genotypes as output from [sim_geno()].
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: logical values, or character
#' string IDs. Numbers are interpreted as character string IDs.
#' @param ... Ignored.
#'
#' @return The input imputed genotypes, with the selected
#' individuals and/or chromsomes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' \dontshow{grav2 <- grav2[1:8,c(1,2)]}
#' dr <- sim_geno(grav2, n_draws=4)
#' # keep just individuals 1:5, chromosome 2
#' drsub <- dr[1:5,2]
#' # keep just chromosome 2
#' drsub2 <- dr[,2]
subset.sim_geno <-
    function(x, ind=NULL, chr=NULL, ...)
    subset.calc_genoprob(x, ind, chr, ...)

#' @export
#' @rdname subset.sim_geno
`[.sim_geno` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)
