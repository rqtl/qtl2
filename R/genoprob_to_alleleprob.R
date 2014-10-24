# genoprob_to_alleleprob
#' Convert genotype probabilities to allele probabilities
#'
#' Reduce genotype probabilities (as calculated by \code{\link{calc_genoprob}}) to allele probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as calculated from \code{\link{calc_genoprob}}.
#' @param quiet IF \code{FALSE}, print messages on progress.
#'
#' @return List of three-dimensional arrays of probabilities, regarding alleles rather than genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#' allele_probs <- genoprob_to_alleleprob(probs)

genoprob_to_alleleprob <-
    function(probs, quiet=TRUE)
{
    is_x_chr <- attr(probs, "is_x_chr")

    for(i in seq(along=probs)) {
        if(!quiet) message(" - Chr ", names(probs)[i])
        attr_i <- attributes(probs[[i]])
        probs[[i]] <- aperm(.genoprob_to_alleleprob(attr(probs, "crosstype"),
                                                    aperm(probs[[i]], c(3, 1, 2)),
                                                    is_x_chr[i]),
                            c(2, 3, 1))
        dimnames(probs[[i]]) <- attr_i$dimnames
    }

    attr(probs, "alleles") <- TRUE
    probs
}
