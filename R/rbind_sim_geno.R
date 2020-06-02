#' Join genotype imputations for different individuals
#'
#' Join multiple genotype imputation objects, as produced by
#' [sim_geno()], for the same set of markers but for different
#' individuals.
#'
#' @param ... Genotype imputation objects as produced by
#' [sim_geno()]. Must have the same set of markers and
#' genotypes.
#'
#' @return An object of class `"sim_geno"`, like the input; see [sim_geno()].
#'
#' @seealso [cbind.sim_geno()], [sim_geno()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' drawsA <- sim_geno(grav2[1:5,], map, error_prob=0.002, n_draws=4)
#' drawsB <- sim_geno(grav2[6:12,], map, error_prob=0.002, n_draws=4)
#' draws <- rbind(drawsA, drawsB)
#'
#' @export
rbind.sim_geno <-
    function(...)
    rbind.calc_genoprob(...)
