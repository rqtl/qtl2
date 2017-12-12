# calc_grid

#' Calculate indicators of which marker/pseudomarker positions are along a fixed grid
#'
#' Construct vectors of logical indicators that indicate which
#' positions correspond to locations along a grid
#'
#' @md
#'
#' @param map A list of numeric vectors; each vector gives marker
#' positions for a single chromosome.
#' @param step Distance between pseudomarkers and markers; if
#' `step=0` no pseudomarkers are inserted.
#' @param off_end Distance beyond terminal markers in which to insert
#' pseudomarkers.
#' @param tol Tolerance for determining whether a pseudomarker would
#' duplicate a marker position.
#'
#' @return A list of logical (TRUE/FALSE) vectors that indicate, for a
#'     marker/pseudomarker map created by
#'     [insert_pseudomarkers()] with `step`>0 and
#'     `stepwidth="fixed"`, which positions correspond to he
#'     locations along the fixed grid.
#'
#' @details The function [insert_pseudomarkers()], with
#'     `stepwidth="fixed"`, will insert a grid of pseudomarkers,
#'     to a marker map. The present function gives a series of
#'     TRUE/FALSE vectors that indicate which positions fall on the
#'     grid. This is for use with [probs_to_grid()], for
#'     reducing genotype probabilities, calculated with
#'     [calc_genoprob()], to just the positions on the grid.
#'     The main value of this is to speed up genome scan computations
#'     in the case of very dense markers, by focusing on just a grid
#'     of positions rather than on all marker locations.
#'
#' @export
#' @keywords utilities
#' @seealso [insert_pseudomarkers()], [probs_to_grid()],
#'     [map_to_grid()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' gmap_w_pmar <- insert_pseudomarkers(iron$gmap, step=1)
#' grid <- calc_grid(iron$gmap, step=1)
calc_grid <-
function(map, step=0, off_end=0, tol=0.01)
{
    # is input not a map but a cross2 object?
    cl <- class(map)
    if(length(cl)==1 && cl=="cross2") {
        map <- map$gmap
        if(is.null(map))
            stop("Input is a cross2 object, but no genetic map found.")
    }
    if(!is.number(step) || step<0)
        stop("step should be a single non-negative number")
    if(!is.number(off_end) || off_end<0)
        stop("off_end should be a single non-negative number")
    if(!is.number(tol) || tol<0)
        stop("tol should be a single non-negative number")

    grid <- map
    chr <- names(grid) <- names(map)
    for(i in seq(along=map)) {
        grid[[i]] <- calc_grid_onechr(map[[i]], step, off_end, tol,
                                      paste0("c", chr[i], ".loc"))
    }

    grid
}



# grid for pseudomarkers
#
# map = vector of marker positions; names = marker names
# step = distance between pseudomarkers
# off_end = amount to go beyond the terminal markers
# tol = tolerance for determining whether a marker hits a pseudomarker
# pmar_stem = leading part of name of pseudomarkers
calc_grid_onechr <-
function(map, step=0, off_end=0, tol=0.01, pmar_stem="loc")
{
    if(!is.number(step) || step < 0)
        stop("step should be a single non-negative number")
    if(!is.number(off_end) || off_end < 0)
        stop("off_end should be a single non-negative number")
    if(!is.number(tol) || tol < 0)
        stop("tol should be a single non-negative number")

    if(step==0) {
        grid <- rep(TRUE, length(map))
        names(grid) <- names(map)
        return(grid)
    }

    if(any(is.na(map))) stop("map values can't be missing")

    if(tol < 0) stop("tol should be >= 0")

    # locations of pseudomarkers
    pmar_map <- seq(min(map)-off_end, max(map)+off_end, by=step)

    # for each pseudomarker, distance to nearest marker
    d <- abs(outer(map, pmar_map, "-"))
    mind <- apply(d, 2, min)

    # omit pseudomarkers that hit a marker
    to_omit <- (mind < tol)
    pmar_map <- pmar_map[!to_omit]
    if(length(pmar_map) == 0) {
        grid <- rep(TRUE, length(map))
        names(grid) <- names(map)
        return(grid)
    }

    # markers that are on pseudomarker grid
    grid <- rep(FALSE, length(map))
    names(grid) <- names(map)
    if(any(to_omit))
        grid[apply(d[,to_omit,drop=FALSE], 2, which.min)] <- TRUE

    if(is.null(names(pmar_map)))
        names(pmar_map) <- create_pseudomarker_names(pmar_map, step, pmar_stem)

    # logical vector indicating grid
    pmar_grid <- rep(TRUE, length(pmar_map))
    names(pmar_grid) <- names(pmar_map)
    grid <- c(grid, pmar_grid)

    # add pseudomarkers to markers
    map <- c(map, pmar_map)

    # sort the map
    o <- order(map)

    grid[o]
}
