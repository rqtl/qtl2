#' Split vector into batches
#'
#' Split a vector into batches, each no longer than `batch_size` and
#' creating at least `n_cores` batches, for use in parallel
#' calculations.
#'
#' @md
#'
#' @param vec A vector to be split into batches
#' @param batch_size Maximum size for each batch
#' @param n_cores Number of compute cores, to be used as a minimum number of batches.
#'
#' @return A list of vectors, each no longer than `batch_size`, and with at least `n_cores` componenets.
#'
#' @export
#' @importFrom parallel detectCores
#' @seealso [batch_cols()]
#' @examples
#' vec_split <- batch_vec(1:304, 50, 8)
#' vec_split2 <- batch_vec(1:304, 50)
batch_vec <-
    function(vec, batch_size=NULL, n_cores=1)
{
    n <- length(vec)

    if(n_cores==0) n_cores <- parallel::detectCores()

    if(is.null(batch_size)) n_batches <- n_cores
    else {
        n_batches <- ceiling(n/batch_size)
        if(n_batches < n_cores) n_batches <- n_cores
    }

    n_per_batch <- rep(floor(n/n_batches), n_batches)
    d <- n - sum(n_per_batch)
    if(d >= 1)
        n_per_batch[seq_len(d)] <- n_per_batch[seq_len(d)]+1

    result <- split(vec, rep(seq_len(n_batches), n_per_batch))
    names(result) <- NULL
    result
}
