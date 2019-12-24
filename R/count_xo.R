#' Count numbers of crossovers
#'
#' Estimate the numbers of crossovers in each individual on each chromosome.
#'
#' @param geno List of matrices of genotypes (output of [maxmarg()] or [viterbi()])
#' or a list of 3d-arrays of genotypes (output of [sim_geno()]).
#' @param quiet If FALSE, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A matrix of crossover counts, individuals x chromosomes, or
#' (if the input was the output of [sim_geno()]) a
#' 3d-array of crossover counts, individuals x chromosomes x
#' imputations.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[1:20,1:2] # subset to first 20 individuals and to chr 1 and 2}
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, map, error_prob=0.002, map_function="c-f")
#' g <- maxmarg(pr)
#' n_xo <- count_xo(g)
#'
#' # imputations
#' imp <- sim_geno(iron, map, error_prob=0.002, map_function="c-f", n_draws=32)
#' n_xo_imp <- count_xo(imp)
#' # sums across chromosomes
#' tot_xo_imp <- apply(n_xo_imp, c(1,3), sum)
#' # mean and SD across imputations
#' summary_xo <- cbind(mean=rowMeans(tot_xo_imp),
#'                     sd=apply(tot_xo_imp, 1, sd))
#'
#' @seealso [locate_xo()]
#'
#' @importFrom stats setNames
#' @export

count_xo <-
    function(geno, quiet=TRUE, cores=1)
{
    if(is.cross2(geno))
        stop('Input geno is a "cross2" object but should be genotypes as from viterbi, maxmarg, or sim_geno')

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    crosstype <- attr(geno, "crosstype")
    if(is.null(crosstype))
        stop("Input geno needs to include a crosstype attribute.")
    is_x_chr <- attr(geno, "is_x_chr")
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(geno))

    # the case of 3d-arrays from sim_geno
    if(length(dim(geno[[1]])) == 3) {
        nind <- vapply(geno, nrow, 1)
        if(length(unique(nind)) > 1)
            stop("Input chromosomes have varying numbers of individuals")
        nchr <- length(geno)
        nind <- nind[1]
        ndraws <- vapply(geno, function(a) dim(a)[3], 1)
        if(length(unique(ndraws)) > 1)
            stop("Input chromosomes have varying numbers of imputations")
        ndraws <- ndraws[1]

        result <- array(dim=stats::setNames(c(nind, ndraws, nchr), NULL))
        dimnames(result) <- list(rownames(geno[[1]]), NULL, names(geno))

        by_chr_func <- function(chr) {
            result <- .count_xo_3d(aperm(geno[[chr]], c(2,1,3)), crosstype, is_x_chr[chr])
            names(result) <- rownames(geno[[chr]])
            result
        }

        result_list <- cluster_lapply(cores, seq(along=geno), by_chr_func)

        for(i in seq(nchr))
            result[,,i] <- result_list[[i]]

        return(aperm(result, c(1,3,2)))
    }

    # the case of matrices from viterbi or maxmarg
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
