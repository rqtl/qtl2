#' Calculate entropy of genotype probability distribution
#'
#' For each individual at each genomic position, calculate the entropy
#' of the genotype probability distribution, as a quantitative summary
#' of the amount of missing information.
#'
#' @md
#'
#' @param probs Genotype probabilities, as calculated from
#' [calc_genoprob()].
#' @param quiet IF `FALSE`, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list of matrices (each matrix is a chromosome and is arranged as individuals x markers).
#'
#' @details We calculate -sum(p log_2 p), where we take 0 log 0 = 0.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' probs <- calc_genoprob(grav2, error_prob=0.002)
#' e <- calc_entropy(probs)
#' e <- do.call("cbind", e) # combine chromosomes into one big matrix
#'
#' # summarize by individual
#' hist(rowMeans(e), breaks=25, main="Ave entropy by individual", xlab="Entropy")
#'
#' # summarize by marker
#' plot(colMeans(e), xlab="marker index", ylab="Average entropy")
calc_entropy <-
    function(probs, quiet=TRUE, cores=1)
{
    if("cross2" %in% class(probs))
        stop('Input probs is a "cross2" object but should be genotype probabilities, as from calc_genoprob')

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # function that does the work
    entropy <-
        function(p, tol=1e-256) {
            p <- p[p>tol]
            -sum(p * log2(p))
        }

    result <- cluster_lapply(cores, seq_along(probs),
                             function(i) apply(probs[[i]], c(1,3), entropy))
    names(result) <- names(probs)
    result
}
