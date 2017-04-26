#' Count numbers of crossovers
#'
#' Estimate the numbers of crossovers in each individual on each chromosome.
#'
#' @param geno List of matrices of genotypes (output of \code{\link{maxmarg}} or \code{\link{viterbi}}).
#' @param quiet If FALSE, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A matrix of crossover counts, individuals x chromosomes
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, map, error_prob=0.002, map_function="c-f")
#' g <- maxmarg(pr)
#' count_xo(g)
#'
#' @export

count_xo <-
    function(geno, quiet=TRUE, cores=1)
{
    if("cross2" %in% class(geno))
        stop('Input geno is a "cross2" object but should be genotypes as from viterbi or maxmarg')

    crosstype <- attr(geno, "crosstype")
    if(is.null(crosstype))
        stop("Input geno needs to include a crosstype attribute.")
    is_x_chr <- attr(geno, "is_x_chr")
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(geno))

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    by_chr_func <- function(chr) {
        result <- .count_xo(t(geno[[chr]]), crosstype, is_x_chr[chr])
        names(result) <- rownames(geno[[chr]])
        result
    }

    result_list <- cluster_lapply(cores, seq(along=geno), by_chr_func)

    result <- do.call("cbind", result_list)
    colnames(result) <- names(geno)
    result
}
