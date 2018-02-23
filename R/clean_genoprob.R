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
#' @param value_threshold Probabilities below this value will be set
#'     to 0.
#' @param column_threshold For genotype columns where the maximum
#'     value is below this threshold, all values will be set to 0.
#'     This must be less than \eqn{1/k} where \eqn{k} is the number of
#'     genotypes.
#' @param ind Optional vector of individuals (logical, numeric, or
#'     character). If provided, only the genotype probabilities for
#'     these individuals will be cleaned, though the full set will be
#'     returned.
#' @param cores Number of CPU cores to use, for parallel calculations.
#'     (If `0`, use [parallel::detectCores()].)
#'     Alternatively, this can be links to a set of cluster sockets, as
#'     produced by [parallel::makeCluster()].
#'
#' @return A cleaned version of the input genotype probabilities
#'     object, `probs`.
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
#' If `ind` is provided, the function is applied only to the
#' designated subset of individuals. This may be useful when only a
#' subset of individuals have been phenotyped, as you may want to zero
#' out genotype columns where that subset of individuals has only
#' negligible probability values.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c("19", "X")] # subset to chr 19 and X}
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, error_prob=0.002)
#'
#' # clean the genotype probabilities
#' # (doesn't really do anything in this case, because there are no small but non-zero values)
#' probs_clean <- clean_genoprob(probs)
#'
#' # clean only the females' genotype probabilities
#' probs_cleanf <- clean_genoprob(probs, ind=names(iron$is_female)[iron$is_female])

#' @export
clean_genoprob <-
    function(probs, value_threshold=1e-6, column_threshold=0.01, ind=NULL, cores=1)
{
    if(!is.null(ind)) {
        # clean the selected subset of individuals
        probs_sub <- clean_genoprob(subset(probs, ind=ind), value_threshold=value_threshold,
                                    column_threshold=column_threshold, cores=cores)

        # paste the result into the main object
        for(i in seq_along(probs))
            probs[[i]][rownames(probs_sub[[i]]),,] <- probs_sub[[i]]

        # exit the function
        return(probs)
    }

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
