#' Calculate heterozygosities
#'
#' Calculate heterozygosites, by individual or by marker
#'
#' @param probs List of arrays of genotype probabilities, as
#'     calculated by [calc_genoprob()].
#' @param by Whether to summarize by individual or marker
#'
#' @param omit_x If TRUE, omit the X chromosome.
#'
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return
#' The result is a vector of estimated heterozygosities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' p <- calc_genoprob(iron, err=0.002)
#'
#' # heterozygosities by individual
#' het_ind <- calc_het(p)
#'
#' # heterozygosities by marker
#' het_mar <- calc_het(p, "marker")
#'
#' @export
calc_het <-
    function(probs, by=c("individual", "marker"), omit_x=TRUE, cores=1)
{
    by <- match.arg(by)
    if(is.null(probs)) stop("probs is NULL")
    if(is.cross2(probs))
        stop('Input probs is a "cross2" object but should be genotype probabilities, as from calc_genoprob')

    aprob_attr <- attr(probs, "alleleprobs")
    if(!is.null(aprob_attr) && aprob_attr) {
        stop("Input probs should not be allele dosages")
    }

    # set up cluster
    cores <- setup_cluster(cores, quiet=TRUE)

    is_x_chr <- attr(probs, "is_x_chr")
    if(any(is_x_chr) && !all(is_x_chr) && omit_x) {
        probs <- probs[,!is_x_chr]
        is_x_chr <- is_x_chr[!is_x_chr]
        omit_x <- FALSE
    }

    n_chr <- length(probs)

    # determine which columns are het
    het_col <- vector("list", n_chr)
    geno <- dimnames(probs)[[2]]
    for(chr in seq_len(n_chr)) {
        a1 <- substr(geno[[chr]], 1, 1)
        a2 <- substr(geno[[chr]], 2, 2)
        if(is_x_chr[chr]) het_col[[chr]] <- (a1 != a2 & a2 != "Y")
        else het_col[[chr]] <- (a1 != a2)
    }

    if(by=="individual") {
        # total markers
        total_mar <- sum(dim(probs)[3,])

        # summarize each chromosome
        result <- cluster_lapply(cores, seq_len(n_chr), function(chr) apply(probs[[chr]][,het_col[[chr]],,drop=FALSE], 1, sum))

        if(length(result)>1) {
            for(i in seq_along(result)[-1])
                result[[1]] <- result[[1]] + result[[i]]
        }
        return(result[[1]]/total_mar)
    }

    # else: by marker
    unlist(cluster_lapply(cores, seq_len(n_chr), function(chr) apply(probs[[chr]][,het_col[[chr]],,drop=FALSE], 3, sum)/
                                                          nrow(probs[[chr]])))
}
