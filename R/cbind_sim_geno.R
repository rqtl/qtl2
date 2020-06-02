#' Join genotype imputations for different chromosomes
#'
#' Join multiple genotype imputation objects, as produced by
#' [sim_geno()], for the same individuals but different
#' chromosomes.
#'
#' @param ... Genotype imputation objects as produced by
#' [sim_geno()]. Must have the same set of individuals.
#'
#' @return An object of class `"sim_geno"`, like the input; see [sim_geno()].
#'
#' @seealso [rbind.sim_geno()], [sim_geno()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' drawsA <- sim_geno(grav2[1:5,1:2], map, error_prob=0.002, n_draws=4)
#' drawsB <- sim_geno(grav2[1:5,3:4], map, error_prob=0.002, n_draws=4)
#' draws <- cbind(drawsA, drawsB)
#'
#' @export
cbind.sim_geno <-
    function(...)
    cbind.calc_genoprob(...)
