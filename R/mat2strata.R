#' Define strata based on rows of a matrix
#'
#' Use the rows of a matrix to define a set of strata for a stratified permutation test
#'
#' @param mat A covariate matrix, as individuals x covariates
#'
#' @return A vector of character strings: for each row of \code{mat},
#' we use \code{\link[base]{paste}} with \code{collapse="|"}.
#'
#' @examples
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' Xcovar <- get_x_covar(iron)
#' perm_strata <- mat2strata(Xcovar)
#'
#' @seealso \code{\link[qtl2geno]{get_x_covar}}, \code{\link{scan1perm}}
#' @export
mat2strata <-
    function(mat)
{
    if(is.null(mat)) return(NULL)

    result <- apply(mat, 1, paste, collapse="|")
    names(result) <- rownames(mat)

    result
}
