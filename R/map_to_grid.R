# map_to_grid
#' Subset a map to positions on a grid
#'
#' Subset a map object to the locations on some grid.
#'
#' @param map A list of vectors of marker positions.
#' @param grid A list of vectors of logical values (aligned with
#' \code{map}), with TRUE indicating the position is on the grid.
#' If missing, we look in \code{map} for an attribute
#' \code{"grid"}.
#'
#' @return Same list as input, but subset to just include
#' pseudomarkers along a grid. The map attribute is similarly subset.
#'
#' @details This is generally for the case of a map created with
#' \code{\link{insert_pseudomarkers}} with \code{step}>0 and
#' \code{stepwidth="fixed"}, so that the pseudomarkers form a grid
#' along each chromosome.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map_w_pmar <- insert_pseudomarkers(grav2$gmap, step=1)
#' sapply(map_w_pmar, length)
#' map_sub <- map_to_grid(map_w_pmar)
#' sapply(map_sub, length)

map_to_grid <-
    function(map, grid=NULL)
{
    if(is.null(grid))
        grid <- attr(map, "grid")

    if(length(map) != length(grid))
        stop("length(grid) [", length(grid), "] != length(map) [", length(map), "]")

    for(i in seq(along=map)) {
        if(is.null(grid[[i]]) || all(grid[[i]])) next

        if(length(map[[i]]) != length(grid[[i]]))
            stop("length(grid) [", length(grid[[i]]), "] != length(map) [",
                 length(map[[i]]), "] for chr ", names(map)[i])

        # subset map
        map[[i]] <- map[[i]][grid[[i]]]
    }

    map
}
