# genoprob_to_alleleprob
#' Convert genotype probabilities to allele probabilities
#'
#' Reduce genotype probabilities (as calculated by
#' \code{\link{calc_genoprob}}) to allele probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as
#' calculated from \code{\link{calc_genoprob}}.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return List of three-dimensional arrays of probabilities,
#' regarding alleles rather than genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#' allele_probs <- genoprob_to_alleleprob(probs)

genoprob_to_alleleprob <-
    function(probs, quiet=TRUE, cores=1)
{
    is_x_chr <- attr(probs, "is_x_chr")

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores.")
        quiet <- TRUE # no more messages
    }

    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        attr_chr <- attributes(probs[[chr]])
        probs[[chr]] <- aperm(.genoprob_to_alleleprob(attr(probs, "crosstype"),
                                                    aperm(probs[[chr]], c(2, 1, 3)), # reorg -> geno x ind x pos
                                                    is_x_chr[chr]),
                            c(2, 1, 3)) # reorg back to ind x geno x pos
        dimnames(probs[[chr]]) <- attr_chr$dimnames
        probs[[chr]]
    }

    chrs <- seq(along=probs)
    probs_attr <- attributes(probs)

    probs <- cluster_lapply(cores, chrs, by_chr_func) # if cores==1, this uses lapply()

    for(at in names(probs_attr))
        attr(probs, at) <- probs_attr[[at]]
    attr(probs, "alleles") <- TRUE
    probs
}
