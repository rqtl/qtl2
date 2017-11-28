#' Batch columns by pattern of missing values
#'
#' Identify batches of columns of a matrix that have the same pattern
#' of missing values.
#'
#' @md
#'
#' @param mat A numeric matrix
#' @param max_batch Maximum batch size
#'
#' @return A list containing the batches, each with two components:
#' `cols` containing numeric indices of the columns in the
#' corresponding batch, and `omit` containing a vector of row indices
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
    function(mat, max_batch=NULL)
{
    if(is.null(max_batch) || max_batch<=0)
        max_batch <- ncol(mat)

    mat <- !is.finite(mat)
    n <- nrow(mat)
    all_true <- rep(TRUE, n)
    result <- NULL

    # first pull out the columns that have 0 NAs
    n_na <- colSums(mat)
    no_na <- (n_na==0)
    if(any(no_na))
        result <- list(list(cols=which(no_na),
                                omit=numeric(0)))

    # now deal with the columns that have exactly 1 NA
    one_na <- (n_na==1)
    if(any(one_na)) {
        wh <- apply(mat[,one_na,drop=FALSE], 2, which)
        spl <- split(which(one_na), wh)
        part2 <- lapply(seq_along(spl), function(i)
                        list(cols=as.numeric(spl[[i]]),
                             omit=as.numeric(names(spl)[i])))

        if(is.null(result)) result <- part2
        else result <- c(result, part2)
    }

    # now the columns with >1 NA
    other_cols <- !(no_na | one_na)
    if(any(other_cols)) {
        other_cols <- seq_len(ncol(mat))[other_cols]

        # pattern of missing data (as character string with 0's and 1's
        pat <- apply(mat[,other_cols,drop=FALSE], 2, function(a) paste(which(a), collapse=":"))

        # the unique patterns
        u <- unique(pat)

        # the result
        part3 <- lapply(u, function(a)
                        list(cols=other_cols[pat==a],
                             omit=as.numeric(strsplit(a, ":")[[1]])))

    } else part3 <- NULL

    if(is.null(result)) result <- part3
    else result <- c(result, part3)

    # reduce to max_batch
    sizes <- vapply(result, function(a) length(a$cols), 1)
    big <- sizes > max_batch
    if(!any(big)) return(result)

    toreduce <- which(big)
    for(i in toreduce) {
        omit <- result[[i]]$omit
        col_list <- batch_vec(result[[i]]$cols, max_batch)
        this <- lapply(col_list, function(a) list(cols=a, omit=omit))
        result <- c(result, this)
    }

   result[-toreduce]
}


# split a vector into batches, each no longer than batch_size
batch_vec <-
    function(vec, batch_size)
{
    n <- length(vec)
    n_batches <- ceiling(n/batch_size)
    n_per_batch <- rep(floor(n/n_batches), n_batches)
    d <- n - sum(n_per_batch)
    if(d >= 1)
        n_per_batch[seq_len(d)] <- n_per_batch[seq_len(d)]+1

    result <- split(vec, rep(seq_len(n_batches), n_per_batch))
    names(result) <- NULL
    result
}
