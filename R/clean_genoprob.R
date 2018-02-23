#' Clean genotype probabilities
#'
#' Clean up genotype probabilities by setting small values to 0 and
#' for a genotype column where the maximum value is rather small, set
#' all values in that column to 0.
#'
#' @md
#'
#' @param probs Genotype probabilities as calculated by
#'     [calc_genoprob()].
#' @param value_threshold Probabilities below this value will be set to 0.
#' @param column_threshold For genotype columns where the maximum
#'     value is below this threshold, all values will be set to 0.
#'     This must be less than \eqn{1/k} where \eqn{k} is the number of genotypes.
#' @param cores Number of CPU cores to use, for parallel calculations.
#'     (If `0`, use [parallel::detectCores()].)
#'     Alternatively, this can be links to a set of cluster sockets, as
#'     produced by [parallel::makeCluster()].
#'
#' @return A cleaned version of the input genotype probabilities object, `probs`.
#'
#' @details
#' In cases where a particular genotype is largely absent,
#' `scan1coef()` and `fit1()` can give unstable estimates of the
#' genotype effects. Cleaning up the genotype probabilities by setting
#' small values to 0 helps to ensure that such effects get set to
#' `NA`.
#'
#' At each position and for each genotype column, we find the maximum
#' probability across individuals. If that maximum is <
#' `column_threshold`, all values in that genotype column at that
#' position are set to 0.
#'
#' In addition, any genotype probabilties that are < `value_threshold`
#' (generally < `column_threshold`) are set to 0.
#'
#' The probabilities are then re-scaled so that the probabilities for
#' each individual at each position sum to 1.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c("19", "X")] # subset to chr 19 and X}
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, error_prob=0.002)
#'
#' # clean the genotype probabilities and paste over original values
#' # (doesn't really do anything in this case, because there are no small but non-zero values)
#' probs <- clean_genoprob(probs)
#'
#' @export
clean_genoprob <-
    function(probs, value_threshold=1e-6, column_threshold=0.01, cores=1)
{
    attrib <- attributes(probs)

    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, seq_along(probs),
                             function(i) {
        this_result <- .clean_genoprob(probs[[i]], value_threshold, column_threshold)
        dimnames(this_result) <- dimnames(probs[[i]])
        this_result })

    for(a in names(attrib)) {
        attr(result, a) <- attrib[[a]]
    }

    result
}
