#' Batch columns by pattern of missing values
#'
#' Identify batches of columns of a matrix that have the same pattern
#' of missing values.
#'
#' @param mat A numeric matrix
#'
#' @return A list containing the batches, each with two components:
#' \code{cols} containing numeric indices of the columns in the
#' corresponding batch, and \code{keep} containing a vector of row indices
#' that have missing values in this batch.
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
    mat <- is.na(mat)
    n <- nrow(mat)
    all_true <- rep(TRUE, n)
    result <- NULL

    # first pull out the columns that have 0 NAs
    n_na <- colSums(mat)
    no_na <- (n_na==0)
    if(any(no_na))
        result <- list(list(cols=which(no_na),
                                keep=numeric(0)))

    # now deal with the columns that have exactly 1 NA
    one_na <- (n_na==1)
    if(any(one_na)) {
        wh <- apply(mat[,one_na,drop=FALSE], 2, which)
        spl <- split(which(one_na), wh)
        part2 <- lapply(seq(along=spl), function(i)
                        list(cols=as.numeric(spl[[i]]),
                             keep=as.numeric(names(spl)[i])))

        if(is.null(result)) result <- part2
        else result <- c(result, part2)
    }

    # now the columns with >1 NA
    other_cols <- !(no_na | one_na)
    if(any(other_cols)) {
        other_cols <- (1:ncol(mat))[other_cols]

        # pattern of missing data (as character string with 0's and 1's
        pat <- apply(mat[,other_cols,drop=FALSE], 2, function(a) paste(which(a), collapse=":"))

        # the unique patterns
        u <- unique(pat)

        # the result
        part3 <- lapply(u, function(a)
                        list(cols=other_cols[pat==a],
                             keep=as.numeric(strsplit(a, ":")[[1]])))

        if(is.null(result)) return(part3)
    } else part3 <- NULL

    if(is.null(part3)) return(result)

    c(result, part3)
}
