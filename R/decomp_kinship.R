#' Calculate eigen decomposition of kinship matrix
#'
#' Calculate the eigen decomposition of a kinship matrix, or of a list of such matrices.
#'
#' @param kinship A square matrix, or a list of square matrices.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return The eigen values and the \bold{transposed} eigen vectors,
#' as a list containing a vector \code{values} and a matrix
#' \code{vectors}.
#'
#' @details The result contains an attribute \code{"eigen_decomp"}.
#'
#' @examples
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
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

    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, kinship, Rcpp_eigen_decomp)
    attr(result, "eigen_decomp") <- TRUE
    result
}
