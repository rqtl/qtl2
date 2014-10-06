# insert_pseudomarkers

#' Insert pseudomarkers into a marker map
#'
#' Insert pseudomarkers into a map of genetic markers, for a single chromosome.
#'
#' @param map Numeric vector of marker positions, for a single
#' chromosome.
#' @param step Distance between pseudomarkers and markers; if
#' \code{step=0} no pseudomarkers are inserted.
#' @param off_end Distance beyond terminal markers in which to insert
#' pseudomarkers.
#' @param stepwidth Indicates whether to use a fixed grid
#' (\code{stepwidth="fixed"}) or to use the maximal distance between
#' pseudomarkers to ensure that no two adjacent markers/pseudomarkers
#' are more than \code{step} apart.
#' @param tol Tolerance for determining whether a pseudomarker would duplicate a marker position.
#' @param pmar_step Character string to serve as the stem for naming the pseudomarkers.
#'
#' @return A vector of positions of pseudomarkers and markers: the
#' input \code{map} vector with pseudomarker positions added. An
#' attribute \code{"index"} is an integer vector that indicates which
#' positions are pseudomarkers (value 0) and which are markers
#' (positive values, indicating the marker indices). If
#' \code{stepwidth=fixed}, a further attributed (\code{"grid"}) is a
#' logical vector that indicates which positions correspond to the
#' fixed grid.
#'
#' @details If \code{stepwidth="fixed"}, a grid of pseudomarkers is added to the marker map.
#'
#' If \code{stepwidth="max"}, a minimal set of pseudomarkers are
#' added, so that the maximum distance between adjacent markers or
#' pseudomarkers is at least \code{step}. If two adjacent markers are
#' separated by less than \code{step}, no pseudomarkers will be added
#' to the interval. If they are more then \code{step} apart, a set of
#' equally-spaced pseudomarkers will be added.
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{calc_genoprob}}
#'
#' @examples
#' library(qtl)
#' data(hyper)
#' chr4map <- pullMap(hyper, chr=4)
#' pmap <- insert_pseudomarkers(pmap, step=1, pmar_stem="c4.loc")
insert_pseudomarkers <-
function(map, step=0, off_end=0, stepwidth=c("fixed", "max"),
         tol=0.01, pmar_stem="loc")
{
   switch(match.arg(stepwidth),
          fixed=insert_pseudomarkers_grid(map, step, off_end, tol, pmar_stem),
          max=insert_pseudomarkers_minimal(map, step, off_end, tol, pmar_stem))
}

# pseudomarker map as grid, ignoring markers
#
# map = vector of marker positions; names = marker names
# step = distance between pseudomarkers
# off_end = amount to go beyond the terminal markers
# tol = tolerance for determining whether a marker hits a pseudomarker
# pmar_stem = leading part of name of pseudomarkers
insert_pseudomarkers_grid <-
function(map, step, off_end=0, tol=0.01, pmar_stem="loc")
{
    if(any(is.na(map))) stop("map values can't be missing")
    if(step==0) {
        attr(map, "index") <- seq(along=map)
        map <- add_pmap_attr(map, "fixed", step, off_end)
        return(map)
    }
    if(step < 0) stop("step should be >= 0")
    if(tol >= step) stop("tol should be << step")
    if(off_end < 0) stop("off_end should be >= 0")
    if(tol < 0) stop("tol should be >= 0")

    # locations of pseudomarkers
    pmar <- seq(min(map)-off_end, max(map)+off_end, by=step)

    # for each pseudomarker, distance to nearest marker
    d <- abs(outer(map, pmar, "-"))
    mind <- apply(d, 2, min)

    # omit pseudomarkers that hit a marker
    to_omit <- (mind < tol)
    pmar <- pmar[!to_omit]
    if(length(pmar) == 0) {
        attr(map, "index") <- seq(along=map)
        attr(map, "grid") <- rep(TRUE, length(map))
        map <- add_pmap_attr(map, "fixed", step, off_end)
        return(map)
    }

    # markers that are on pseudomarker grid
    grid <- rep(FALSE, length(map))
    if(any(to_omit))
        grid[apply(d[,to_omit,drop=FALSE], 2, which.min)] <- TRUE

    # no. digits to use in name
    digits <- ceiling(-log10(step))
    digits <- ifelse(digits < 0, 0, digits)

    names(pmar) <- paste0(pmar_stem, round(pmar, digits))

    # index of markers vs pseudomarkers
    index <- c(seq(along=map), rep(0, length(pmar)))
    # logical vector indicating grid
    grid <- c(grid, rep(TRUE, length(pmar)))

    # add pseudomarkers to markers
    map <- c(map, pmar)

    # sort the map
    o <- order(map)
    map <- map[o]

    attr(map, "index") <- index[o]
    attr(map, "grid") <- grid[o]
    map <- add_pmap_attr(map, "fixed", step, off_end)

    map
}

# pseudomarker map, minimal number of pseudomarkers to add
#
# map = vector of marker positions; names = marker names
# step = min distance between markers or pseudomarkers
# off_end = amount to go beyond the terminal markers
# tol = tolerance for determining whether a marker hits a pseudomarker
# pmar_stem = leading part of name of pseudomarkers
insert_pseudomarkers_minimal <-
function(map, step, off_end=0, tol=0.01, pmar_stem="loc")
{
    if(any(is.na(map))) stop("map values can't be missing")
    if(step==0) {
        attr(map, "index") <- seq(along=map)
        map <- add_pmap_attr(map, "max", step, off_end)
        return(map)
    }
    if(step < 0) stop("step should be >= 0")
    if(tol >= step) stop("tol should be << step")
    if(off_end < 0) stop("off_end should be >= 0")
    if(tol < 0) stop("tol should be >= 0")

    # distance between markers
    d <- diff(map)
    if(any(d < 0))
        stop("marker positions should be sorted")

    pmar <- c(seq(min(map), min(map)-off_end, by=-step)[-1],
              seq(max(map), max(map)+off_end, by=step)[-1])

    if(length(map) == 1) { # just one marker
        index <- c(1, rep(0, length(pmar)))
        map <- c(map, pmar)
        o <- order(map)
        map <- map[o]
        attr(map, "index") <- index[o]
        map <- add_pmap_attr(map, "max", step, off_end)
        return(map)
    }

    # intervals needing pseudomarkers (d = distances between adjacent markers
    int2fill <- which(d > step)

    pmar <- c(pmar, unlist(lapply(int2fill, function(a) {
        n <- ceiling((map[a+1]-map[a])/step)+1
        seq(map[a], map[a+1], length.out=n)[-c(1, n)] })))

    if(length(pmar) == 0) {
        attr(map, "index") <- seq(along=map)
        map <- add_pmap_attr(map, "max", step, off_end)
        return(map)
    }

    # no. digits to use in name
    digits <- ceiling(-log10(step))
    digits <- ifelse(digits < 0, 0, digits)

    names(pmar) <- paste0(pmar_stem, round(pmar, digits))

    # index of markers vs pseudomarkers
    index <- c(seq(along=map), rep(0, length(pmar)))

    # add pseudomarkers to markers
    map <- c(map, pmar)

    # sort the map
    o <- order(map)
    map <- map[o]
    map <- add_pmap_attr(map, "max", step, off_end)
    attr(map, "index") <- index[o]

    map
}

add_pmap_attr <-
function(map, stepwidth, step, off_end)
{
    attr(map, "stepwidth") <- stepwidth
    attr(map, "step") <- step
    attr(map, "off_end") <- off_end
    map
}
