# insert_pseudomarkers

#' Insert pseudomarkers into a marker map
#'
#' Insert pseudomarkers into a map of genetic markers, for a single chromosome.
#'
#' @md
#'
#' @param map A list of numeric vectors; each vector gives marker
#' positions for a single chromosome.
#' @param step Distance between pseudomarkers and markers; if
#' `step=0` no pseudomarkers are inserted.
#' @param off_end Distance beyond terminal markers in which to insert
#' pseudomarkers.
#' @param stepwidth Indicates whether to use a fixed grid
#' (`stepwidth="fixed"`) or to use the maximal distance between
#' pseudomarkers to ensure that no two adjacent markers/pseudomarkers
#' are more than `step` apart.
#' @param pseudomarker_map A map of pseudomarker locations; if provided the
#' `step`, `off_end`, and `stepwidth` arguments are
#' ignored.
#' @param tol Tolerance for determining whether a pseudomarker would duplicate a marker position.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list like the input `map` with pseudomarkers
#' inserted. Will also have an attribute `"is_x_chr"`, taken
#' from the input `map`.
#'
#' @details If `stepwidth="fixed"`, a grid of pseudomarkers is
#' added to the marker map.
#'
#' If `stepwidth="max"`, a minimal set of pseudomarkers are
#' added, so that the maximum distance between adjacent markers or
#' pseudomarkers is at least `step`. If two adjacent markers are
#' separated by less than `step`, no pseudomarkers will be added
#' to the interval. If they are more then `step` apart, a set of
#' equally-spaced pseudomarkers will be added.
#'
#' If `pseudomarker_map` is provided, then the `step`,
#' `off_end`, and `stepwidth` arguments are ignored, and the
#' input `pseudomarker_map` is taken to be the set of
#' pseudomarker positions.
#'
#' @export
#' @keywords utilities
#' @seealso [calc_genoprob()], [calc_grid()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' gmap_w_pmar <- insert_pseudomarkers(iron$gmap, step=1)
insert_pseudomarkers <-
function(map, step=0, off_end=0, stepwidth=c("fixed", "max"),
         pseudomarker_map=NULL, tol=0.01, cores=1)
{
    if(is.null(map))
        stop("Input map is NULL")

    # is input not a map but a cross2 object?
    cl <- class(map)
    if(length(cl)==1 && cl=="cross2") {
        map <- map$gmap
        if(is.null(map))
            stop("Input is a cross2 object, but no genetic map found.")
    }

    # check inputs
    stepwidth <- match.arg(stepwidth)
    if(is.null(pseudomarker_map)) {
        pseudomarker_map <- vector("list", length(map))
    }
    else {
        if(!is.list(pseudomarker_map))
            stop("pseudomarker_map should be a list")
        if(length(pseudomarker_map) != length(map))
            stop("length(pseudomarker_map) != length(map)")
    }
    if(!is.number(step) || step<0)
        stop("step should be a single non-negative number")
    if(!is.number(off_end) || off_end<0)
        stop("off_end should be a single non-negative number")
    if(!is.number(tol) || tol<0)
        stop("tol should be a single non-negative number")

    chr <- names(map)

    # set up cluster
    cores <- setup_cluster(cores)

    # do the work, possible with multiple CPU
    newmap <- cluster_lapply(cores, seq(along=map), function(i)
        insert_pseudomarkers_onechr(map[[i]], step=step, off_end=off_end, stepwidth=stepwidth,
                                    pseudomarker_map[[i]], tol=tol, paste0("c", chr[i], ".loc")))
    names(newmap) <- chr

    attr(newmap, "is_x_chr") <- attr(map, "is_x_chr")
    newmap
}



insert_pseudomarkers_onechr <-
function(map, step=0, off_end=0, stepwidth=c("fixed", "max"),
         pseudomarker_map=NULL, tol=0.01, pmar_stem="loc")
{
    if(any(is.na(map))) stop("map values can't be missing")

    if(tol < 0) stop("tol should be >= 0")

    if(is.null(pseudomarker_map)) {
        stepwidth <- match.arg(stepwidth)
        if(step < 0) stop("step should be >= 0")
        if(step > 0 && tol >= step) stop("tol should be << step")
        if(off_end < 0) stop("off_end should be >= 0")
    }
    else {
        if(any(is.na(pseudomarker_map)))
            stop("pseudomarker_map values can't be missing")
        stepwidth <- "custom"
    }

    switch(stepwidth,
           fixed=insert_pseudomarkers_grid(map, step, off_end, tol, pmar_stem),
           max=insert_pseudomarkers_minimal(map, step, off_end, tol, pmar_stem),
           custom=combine_markers_with_pseudomarkers(map, pseudomarker_map, tol, pmar_stem))
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
    if(step==0) return(map)

    # locations of pseudomarkers
    pmar <- seq(min(map)-off_end, max(map)+off_end, by=step)

    combine_markers_with_pseudomarkers(map, pmar, tol, pmar_stem, step)
}

# combine the markers in "map" with the pseudomarker positions in
# "pmar". treat the positions in "pmar" as a grid, and return a
# logical vector of which values in output are on the grid.
combine_markers_with_pseudomarkers <-
function(map, pmar_map, tol=0.01, pmar_stem="loc", step=0.1)
{
    # for each pseudomarker, distance to nearest marker
    d <- abs(outer(map, pmar_map, "-"))
    mind <- apply(d, 2, min)

    # omit pseudomarkers that hit a marker
    to_omit <- (mind < tol)
    pmar_map <- pmar_map[!to_omit]
    if(length(pmar_map) == 0) return(map)

    if(is.null(names(pmar_map)))
        names(pmar_map) <- create_pseudomarker_names(pmar_map, step, pmar_stem)

    # add pseudomarkers to markers and sort
    sort(c(map, pmar_map))
}

# create names for pseudomarkers
create_pseudomarker_names <-
function(pmar_map, step, pmar_stem="loc")
{
    digits <- ceiling(-log10(step))
    digits <- ifelse(digits < 0, 0, digits)

    pmar_nam <- paste0(pmar_stem, round(pmar_map, digits))

   while(any_duplicates(pmar_nam)) {
       digits <- digits + 1
       pmar_nam <- paste0(pmar_stem, round(pmar_map, digits))
   }
   pmar_nam
}

any_duplicates <-
    function(x)
    length(unique(x)) < length(x)


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
    if(step==0) return(map)

    # distance between markers
    d <- diff(map)
    if(any(d < 0))
        stop("marker positions should be sorted")

    pmar_map <- c(seq(min(map), min(map)-off_end, by=-step)[-1],
                  seq(max(map), max(map)+off_end, by=step)[-1])

    if(length(map) == 1) { # just one marker
        map <- c(map, pmar_map)
        o <- order(map)
        map <- map[o]
        return(map)
    }

    # intervals needing pseudomarkers (d = distances between adjacent markers
    int2fill <- which(d > step)

    pmar_map <- c(pmar_map, unlist(lapply(int2fill, function(a) {
        n <- ceiling((map[a+1]-map[a])/step)+1
        seq(map[a], map[a+1], length.out=n)[-c(1, n)] })))

    if(length(pmar_map) == 0) {
        return(map)
    }

    names(pmar_map) <- create_pseudomarker_names(pmar_map, step, pmar_stem)

    # add pseudomarkers to markers and sort
    sort(c(map, pmar_map))
}
