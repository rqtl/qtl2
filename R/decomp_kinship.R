#' Calculate eigen decomposition of kinship matrix
#'
#' Calculate the eigen decomposition of a kinship matrix, or of a list of such matrices.
#'
#' @md
#'
#' @param kinship A square matrix, or a list of square matrices.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return The eigen values and the **transposed** eigen vectors,
#' as a list containing a vector `values` and a matrix
#' `vectors`.
#'
#' @details The result contains an attribute `"eigen_decomp"`.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#' K <- calc_kinship(probs)
#'
#' Ke <- decomp_kinship(K)
#'
#' @export
decomp_kinship <-
    function(kinship, cores=1)
{
    # already done?
    if(is_kinship_decomposed(kinship))
        return(kinship) # no need to do it again

    if(is.matrix(kinship)) {
        if(ncol(kinship) != nrow(kinship))
            stop("matrix must be square")
        if(ncol(kinship) == 0)
            stop("matrix has dimension (0,0)")
        return(Rcpp_eigen_decomp(kinship))
    }

    if(!is.list(kinship))
        stop("kinship should be either a square matrix or a list of square matrices")

    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, kinship, Rcpp_eigen_decomp)
    attr(result, "eigen_decomp") <- TRUE
    result
}
