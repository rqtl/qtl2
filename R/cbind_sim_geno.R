#' Join genotype imputations for different chromosomes
#'
#' Join multiple genotype imputation objects, as produced by
#' \code{\link{sim_geno}} for the same individuals but different
#' chromosomes.
#'
#' @param ... Genotype imputation objects as produced by
#' \code{\link{sim_geno}}. Must have the same set of individuals.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' drawsA <- sim_geno(grav2[1:5,1:2], step=1, error_prob=0.002, n_draws=4)
#' drawsB <- sim_geno(grav2[1:5,3:4], step=1, error_prob=0.002, n_draws=4)
#' draws <- cbind(drawsA, drawsB)
#'
#' @export
cbind.sim_geno <-
    function(...)
    cbind.calc_genoprob(...)
