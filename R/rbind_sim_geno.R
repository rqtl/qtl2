#' Join genotype imputations for different individuals
#'
#' Join multiple genotype imputation objects, as produced by \code{\link{sim_geno}} for different individuals.
#'
#' @param ... Genotype imputation objects as produced by
#' \code{\link{sim_geno}}. Must have the same set of markers and
#' genotypes.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' drawsA <- sim_geno(grav2[1:5,], step=1, error_prob=0.002, n_draws=4)
#' drawsB <- sim_geno(grav2[6:12,], step=1, error_prob=0.002, n_draws=4)
#' draws <- rbind(drawsA, drawsB)
#'
#' @export
rbind.sim_geno <-
    function(...)
    rbind.calc_genoprob(...)
