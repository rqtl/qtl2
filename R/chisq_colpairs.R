#' Chi-square test on all pairs of columns
#'
#' Perform a chi-square test for independence for all pairs of columns of a matrix.
#'
#' @md
#'
#' @param x A matrix of positive integers. `NA`s and values <= 0 are treated as missing.
#'
#' @return A matrix of size p x p, where p is the number of columns in
#'     the input matrix `x`, containing the chi-square test
#'     statistics for independence, applied to pairs of columns of
#'     `x`. The diagonal of the result will be all `NA`s.
#'
#' @keywords htest
#'
#' @export
#'
#' @examples
#' z <- matrix(sample(1:2, 500, replace=TRUE), ncol=5)
#' chisq_colpairs(z)

chisq_colpairs <-
    function(x)
{
    if(!is.matrix(x))
        x <- as.matrix(x)
    if(ncol(x) < 2)
        stop("ncol(x) should be >= 2")

    result <- .chisq_colpairs(x)
    dimnames(result) <- list(colnames(x), colnames(x))
    diag(result) <- NA

    result
}
