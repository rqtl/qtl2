#' Identify the largest subset of markers that are some distance apart
#'
#' Identify the largest subset of markers for which no two adjacent
#' markers are separated by less than some specified distance; if
#' weights are provided, find the marker subset for which the sum of
#' the weights is maximized.
#'
#' @param map Either a vector of marker positions, or a list of such
#' vectors (one vector per chromosome)
#'
#' @param min_d Minimum distance between markers
#'
#' @param weights An object of the same shape as \code{map}: either a
#' vector of the same length, or a list with the same length and whose
#' components are the same lengths.
#'
#' @return The selected subset of marker positions, either as a vector
#' or a list of vectors, according to the nature of \code{map}.
#'
#' @details Let \eqn{d_i}{d[i]} be the location of marker \eqn{i}, for
#' \eqn{i \in 1, \dots, M}{i in 1, \dots, M}.  We use the dynamic
#' programming algorithm of Broman and Weber (1999) to identify the
#' subset of markers \eqn{i_1, \dots, i_k}{i[1], \dots, i[k]} for
#' which \eqn{d_{i_{j+1}} - d_{i_j} \le}{d(i[j+1]) - d(i[j]) <=}
#' \code{min.distance} and \eqn{\sum w_{i_j}}{sum w(i[j])} is
#' maximized.
#'
#' If there are multiple optimal subsets, we pick one at random.
#'
#' @references Broman KW, Weber JL (1999) Method for constructing
#' confidently ordered linkage maps. Genet Epidemiol 16:337--343.
#'
#' @export
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#'
#' # grap genetic map
#' gmap <- grav2$gmap
#'
#' # subset to markers that are >= 1 cM apart
#' gmap_sub <- pick_marker_subset(gmap, 1)
pick_marker_subset <-
    function(map, min_d=1, weights=NULL)
{
    # multiple chromosomes
    if(is.list(map)) {
        if(!is.null(weights)) {
            if(!is.list(weights) || length(weights) != length(map) ||
               any(vapply(weights, length, 1) != vapply(map, length, 1)))
                stop("map and weights are different shapes")
        }
        for(i in seq_along(map)) {
            if(is.null(weights)) wts <- rep(1, length(map[[i]]))
            else wts <- weights[[i]]

            map[[i]] <- pick_marker_subset(map[[i]], min_d, wts)
        }
        return(map)
    }

    # one chromosome
    if(is.null(weights))
        weights <- rep(1, length(map))
    else if (length(weights) != length(map))
        stop("length(weights) != length(map)")

    # call c++ function
    index <- .pick_marker_subset(map, min_d, weights)

    map[index]
}
