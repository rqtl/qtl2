#' Reduce the lengths of gaps in a map
#'
#' Reduce the lengths of gaps in a map
#'
#' @md
#'
#' @param map Genetic map as a list of vectors (each vector is a
#' chromosome and contains the marker positions).
#' @param min_gap Minimum gap length to return.
#'
#' @return
#' Input map with any gaps greater than `min_gap` reduced to `min_gap`.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' rev_map <- reduce_map_gaps(iron$gmap, 30)
#'
#' @seealso [find_map_gaps()]
#'
#' @export
reduce_map_gaps <-
    function(map, min_gap=50)
{
    for(i in seq_along(map)) {

        # inter-marker distances
        d <- diff(map[[i]])
        d[d > min_gap] <- min_gap

        # back to original positions
        m <- cumsum(c(min(map[[i]]), d))

        # substitute back into object
        names(m) <- names(map[[i]])
        map[[i]] <- m
    }

    map
}
