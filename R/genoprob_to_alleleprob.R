# genoprob_to_alleleprob
#' Convert genotype probabilities to allele probabilities
#'
#' Reduce genotype probabilities (as calculated by
#' \code{\link{calc_genoprob}}) to allele probabilities.
#'
#' @param probs Genotype probabilities, as calculated from
#' \code{\link{calc_genoprob}}.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return The \code{probs} input with probabilities
#' collapsed to alleles rather than genotypes.
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
    # already converted?
    ap <- probs$alleleprobs
    if(!is.null(ap) && ap) return(probs)

    is_x_chr <- probs$is_x_chr

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores.")
        quiet <- TRUE # no more messages
    }

    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", probs$chrID[chr])
        result <- aperm(.genoprob_to_alleleprob(probs$crosstype,
                                                aperm(probs$probs[[chr]], c(2, 1, 3)), # reorg -> geno x ind x pos
                                                is_x_chr[chr]),
                        c(2, 1, 3)) # reorg back to ind x geno x pos

        # allele names
        dn <- dimnames(probs$probs[[chr]])
        dn[[2]] <- probs$alleles
        dimnames(result) <- dn

        result
    }

    chrs <- seq(along=probs$chrID)

    probs$probs <- cluster_lapply(cores, chrs, by_chr_func) # if cores==1, this uses lapply()
    names(probs$probs) <- probs$chrID

    probs$alleleprobs <- TRUE
    probs
}
