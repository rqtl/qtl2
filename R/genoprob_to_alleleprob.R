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

    if("cluster" %in% class(cores) && "SOCKcluster" %in% class(cores)) { # cluster already set
        cluster_ready <- TRUE
        if(!quiet) message(" - Using ", length(cores), " cores.")
        quiet <- TRUE # no more messages
    } else {
        cluster_ready <- FALSE
        if(cores==0) cores <- parallel::detectCores() # if 0, detect cores
        if(cores > 1) {
            if(!quiet) message(" - Using ", cores, " cores.")
            quiet <- TRUE # no more messages
        }
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
    if(!cluster_ready && cores<=1) { # no parallel processing
        probs <- lapply(chrs, by_chr_func)
    }
    else if(cluster_ready || Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        if(!cluster_ready) {
            cores <- parallel::makeCluster(cores)
            on.exit(parallel::stopCluster(cores))
        }
        probs <- parallel::clusterApply(cores, chrs, by_chr_func)
    }
    else {
        probs <- parallel::mclapply(chrs, by_chr_func, mc.cores=cores)
    }

    for(at in names(probs_attr))
        attr(probs, at) <- probs_attr[[at]]
    attr(probs, "alleles") <- TRUE
    probs
}
