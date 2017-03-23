# maxmarg
#' Find genotypes with maximum marginal probabilities
#'
#' For each individual at each position, find the genotype with the maximum marginal probability.
#'
#' @param probs Genotype probabilities, as calculated from
#' \code{\link{calc_genoprob}}.
#' @param map Map of pseudomarkers in \code{probs}. Used only if \code{chr} and \code{pos} are provided.
#' @param minprob Minimum probability for making a call. If maximum
#' probability is less then this value, give \code{NA}.
#' @param chr If provided (along with \code{pos}), consider only the single specified position.
#' @param pos If provided (along with \code{chr}), consider only the single specified position.
#' @param return_char If TRUE, return genotype names as character strings.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return If \code{chr} and \code{pos} are provided, a vector of
#' genotypes is returned. In this case, \code{map} is needed.
#'
#' Otherwise, the result is a object like that returned by \code{\link{viterbi}},
#' A list of two-dimensional arrays of imputed genotypes,
#' individuals x positions. Also includes these attributes:
#' \itemize{
#' \item \code{crosstype} - The cross type of the input \code{cross}.
#' \item \code{is_x_chr} - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input \code{cross}.
#' \item \code{alleles} - Vector of allele codes, from input
#'     \code{cross}.
#' }
#'
#' @seealso \code{\link{sim_geno}}, \code{\link{viterbi}}
#'
#' @export
#'
#' @examples
#' # load data and calculate genotype probabilities
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' pr <- calc_genoprob(iron, error_prob=0.002)
#'
#' # full set of imputed genotypes
#' ginf <- maxmarg(pr)
#'
#' # imputed genotypes at a fixed position
#' g <- maxmarg(pr, iron$gmap, chr=8, pos=45.5)
#'
#' # return genotype names rather than integers
#' g <- maxmarg(pr, iron$gmap, chr=8, pos=45.5, return_char=TRUE)
maxmarg <-
    function(probs, map=NULL, minprob=0.95, chr=NULL, pos=NULL,
             return_char=FALSE, quiet=TRUE, cores=1)
{
    if(!is.null(chr) || !is.null(pos)) {
        if(is.null(chr) || is.null(pos) || is.null(map))
            stop("Provide all of chr, pos, and map, or none of them")
        if(length(chr) > 1 || length(pos) > 1) {
            warning("chr and pos should have length 1; using the first values")
            chr <- chr[1]
            pos <- pos[1]
        }

        # possibly subset the map
        if(length(map) != length(probs) || !all(names(map) == names(probs))) {
            chr <- names(probs)
            if(!all(chr %in% names(map)))
                stop("map doesn't contain all of the necessary chromosomes")
            map <- map[chr]
        }

        chr <- as.character(chr)
        marker <- find_marker(map, chr, pos)

        if(class(probs)[1]=="feather_calc_genoprob")
            probs <- list("1"=probs[,,marker][[1]])
        else
            probs <- list("1"=probs[[chr]][,,marker,drop=FALSE])
        names(probs) <- chr
        class(probs) <- c("calc_genoprob", "list")
    }

    if(return_char && is.null(colnames(probs[[1]]))) {
        warning("No genotype names included in probs; returning integers")
        return_char <- FALSE
    }

    dn <- dimnames(probs)

    # function that does the work
    by_chr_func <- function(chr, return_char) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        result <- .maxmarg(aperm(probs[[chr]], c(2,3,1)), minprob=minprob) # convert to pos x gen x ind
        if(return_char && !is.null(dn[[2]][[chr]]))
            result[,1:ncol(result)] <- dn[[2]][[chr]][result]
        dimnames(result) <- list(dn[[1]], dn[[3]][[chr]])

        result
    }

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # run and combine results
    result <- cluster_lapply(cores, seq_len(length(probs)), by_chr_func, return_char=return_char)
    names(result) <- names(probs)

    # if chr and pos given, return a vector
    if(!is.null(pos)) return(result[[1]][,1])

    # otherwise, return an object like vitebri()
    probs <- result
    class(probs) <- c("viterbi", "list")
    probs
}
