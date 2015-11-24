#' Batch columns by pattern of missing values
#'
#' Identify batches of columns of a matrix that have the same pattern
#' of missing values.
#'
#' @param mat A numeric matrix
#'
#' @return A list containing the batches, each with two components:
#' \code{cols} containing numeric indices of the columns in the
#' corresponding batch, and \code{keep} containing a logical vector
#' with \code{TRUE} indicating that the corresponding row is not
#' missing in this batch.
#'
#' @export
#' @examples
#' x <- rbind(c( 1,  2,  3, 13, 16),
#'            c( 4,  5,  6, 14, 17),
#'            c( 7, NA,  8, NA, 18),
#'            c(NA, NA, NA, NA, 19),
#'            c(10, 11, 12, 15, 20))
#' batch_cols(x)
batch_cols <-
    function(mat)
{
    # pattern of missing data (as character string with 0's and 1's
    pat <- apply(is.na(mat), 2, function(a) paste(as.numeric(a), collapse=""))

    # the unique patterns
    u <- unique(pat)

    # the result
    lapply(u, function(a)
           list(cols=which(pat==a),
                keep=strsplit(a, "")[[1]]=="0"))
}
