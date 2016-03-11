#' Reduce markers to a subset of more-evenly-spaced ones
#'
#' Find the largest subset of markers such that no two adjacent
#' markers are separated by less than some distance.
#'
#' @param map A list with each component being a vector with the
#' marker positions for a chromosome.
#' @param weights A (optional) list of weights on the markers; same
#' size as \code{map}.
#' @param min_distance Minimum distance between markers.
#'
#' @return A list like the input \code{map}, but with the selected
#' subset of markers.
#'
#' @details Uses a dynamic programming algorithm to find, for each
#' chromosome, the subset of markers for with max(\code{weights}) is
#' maximal, subject to the constraint that no two adjacent markers may
#' be separated by more than \code{min_distance}.
#'
#' @references Broman KW, Weber JL (1999) Method for constructing
#' confidently ordered linkage maps. Genet Epidemiol 16:337--343
#'
#' @examples
#' library(qtl)
#'
#' @export
reduce_markers <-
    function(map, weights=NULL, min_distance=1)
{
    if(is.null(weights))
        weights <- lapply(map, function(a) rep(1, length(a)))

    stopifnot(length(map) == length(weights))
    nmar <- vapply(map, length, 1)
    nwts <- vapply(weights, length, 1)
    if(!all(nmar == nwts))
        stop("Different numbers of markers and weights on chr",
             paste(names(map)[nmar != nwts], collapse=" "))

    stopifnot(min_distance > 0)

    for(i in seq(along=map)) {
        if(length(map[[i]]) < 2) next
        map[[i]] <- map[[i]][.reduce_markers(map[[i]], weights[[i]],
                                             min_distance)]
    }

    map
}
