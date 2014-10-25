# genoprob_to_alleleprob
#' Convert genotype probabilities to allele probabilities
#'
#' Reduce genotype probabilities (as calculated by
#' \code{\link{calc_genoprob}}) to allele probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as
#' calculated from \code{\link{calc_genoprob}}.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param n_cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#'
#' @return List of three-dimensional arrays of probabilities,
#' regarding alleles rather than genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#' allele_probs <- genoprob_to_alleleprob(probs)

genoprob_to_alleleprob <-
    function(probs, quiet=TRUE, n_cores=1)
{
    is_x_chr <- attr(probs, "is_x_chr")

    if(n_cores==0) n_cores <- parallel::detectCores() # if 0, detect cores
    if(n_cores > 1) {
        if(!quiet) message(" - Using ", n_cores, " cores.")
        quiet <- TRUE # no more messages
    }

    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        attr_chr <- attributes(probs[[chr]])
        probs[[chr]] <- aperm(.genoprob_to_alleleprob(attr(probs, "crosstype"),
                                                    aperm(probs[[chr]], c(3, 1, 2)),
                                                    is_x_chr[chr]),
                            c(2, 3, 1))
        dimnames(probs[[chr]]) <- attr_chr$dimnames
        probs[[chr]]
    }

    chrs <- seq(along=probs)
    probs_attr <- attributes(probs)
    if(n_cores<=1) { # no parallel processing
        probs <- lapply(chrs, by_chr_func)
    }
    else if(Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))
        probs <- parallel::clusterApply(cl, chrs, by_chr_func)
    }
    else {
        probs <- parallel::mclapply(chrs, by_chr_func, mc.cores=n_cores)
    }

    for(at in names(probs_attr))
        attr(probs, at) <- probs_attr[[at]]
    attr(probs, "alleles") <- TRUE
    probs
}
