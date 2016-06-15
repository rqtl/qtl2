# maxmarg
#' Find genotypes with maximum marginal probabilities
#'
#' For each individual at each position, find the genotype with the maximum marginal probability.
#'
#' @param probs Genotype probabilities, as calculated from
#' \code{\link{calc_genoprob}}.
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
#' @return
#' If \code{chr} and \code{pos} are provided, a vector of genotypes is returned
#'
#' Otherwise, the result is a object like that returned by \code{\link{viterbi}},
#' A list containing the following:
#' \itemize{
#' \item \code{geno} - List of two-dimensional arrays of imputed genotypes,
#' individuals x positions.
#' \item \code{map} - The genetic map as a list of vectors of marker positions.
#' \item \code{grid} - A list of logical vectors, indicating which
#'     positions correspond to a grid of markers/pseudomarkers. (may be
#'     absent)
#' \item \code{crosstype} - The cross type of the input \code{cross}.
#' \item \code{is_x_chr} - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input \code{cross}.
#' \item \code{is_female} - Vector of indicators of which individuals are female, from input
#'     \code{cross}.
#' \item \code{cross_info} - Matrix of cross information for the
#'     individuals, from input \code{cross}.
#' \item \code{alleles} - Vector of allele codes, from input
#'     \code{cross}.
#' \item \code{alleleprob} - Logical value (\code{FALSE}) that
#'     indicates whether the probabilities are compressed to allele
#'     probabilities, as from \code{\link{genoprob_to_alleleprob}}.
#' \item \code{step} - the value of the \code{step} argument.
#' \item \code{off_end} - the value of the \code{off_end} argument.
#' \item \code{stepwidth} - the value of the \code{stepwidth} argument.
#' \item \code{error_prob} - the value of the \code{error_prob} argument.
#' \item \code{map_function} - the value of the \code{map_function} argument.
#' }
#'
#' @seealso \code{\link{sim_geno}}, \code{\link{viterbi}}
#'
#' @export
#'
#' @examples
#' # load data and calculate genotype probabilities
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' pr <- calc_genoprob(iron, step=0, error_prob=0.002)
#'
#' # full set of imputed genotypes
#' ginf <- maxmarg(pr)
#'
#' # imputed genotypes at a fixed position
#' g <- maxmarg(pr, chr=8, pos=45.5)
#'
#' # return genotype names rather than integers
#' g <- maxmarg(pr, chr=8, pos=45.5, return_char=TRUE)
maxmarg <-
    function(probs, minprob=0.95, chr=NULL, pos=NULL,
             return_char=FALSE, quiet=TRUE, cores=1)
{
    if(!is.null(chr) || !is.null(pos)) {
        if(is.null(chr) || is.null(pos))
            stop("Provide both chr and pos, or neither")
        if(length(chr) > 1 || length(pos) > 1) {
            warning("chr and pos should have length 1; using the first values")
            chr <- chr[1]
            pos <- pos[1]
        }

        chr <- as.character(chr)
        marker <- find_marker(probs$map, chr, pos)

        probs$probs <- list("1"=probs$probs[[chr]][,,marker,drop=FALSE])
    }

    if(return_char && is.null(colnames(probs$probs[[1]]))) {
        warning("No genotype names included in probs; returning integers")
        return_char <- FALSE
    }

    # function that does the work
    by_chr_func <- function(chr, return_char) {
        if(!quiet) message(" - Chr ", names(probs$probs)[chr])
        dn <- dimnames(probs$probs[[chr]])
        result <- .maxmarg(aperm(probs$probs[[chr]], c(2,3,1)), minprob=minprob) # convert to pos x gen x ind
        if(return_char && !is.null(dn[[2]]))
            result[,1:ncol(result)] <- dn[[2]][result]
        dimnames(result) <- list(dn[[1]], dn[[3]])
        result
    }

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # run and combine results
    chrs <- seq(along=probs$probs)
    if(n_cores(cores) == 1) {
        result <- lapply(chrs, by_chr_func, return_char=return_char)
    }
    else {
        result <- cluster_lapply(cores, chrs, by_chr_func, return_char=return_char)
    }
    names(result) <- names(probs$probs)

    # if chr and pos given, return a vector
    if(!is.null(pos)) return(result[[1]][,1])

    # otherwise, return an object like vitebri()
    probs$probs <- result
    names(probs)[names(probs)=="probs"] <- "geno"
    class(probs) <- c("viterbi", "list")
    probs
}
