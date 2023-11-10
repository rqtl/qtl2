#' Reduce markers to a subset of more-evenly-spaced ones
#'
#' Find the largest subset of markers such that no two adjacent
#' markers are separated by less than some distance.
#'
#' @param map A list with each component being a vector with the
#' marker positions for a chromosome.
#' @param min_distance Minimum distance between markers.
#' @param weights A (optional) list of weights on the markers; same
#' size as `map`.
#' @param max_batch Maximum number of markers to consider in a batch
#' @param batch_distance_mult If working with batches of markers,
#' reduce `min_distance` by this multiple.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list like the input `map`, but with the selected
#' subset of markers.
#'
#' @details Uses a dynamic programming algorithm to find, for each
#' chromosome, the subset of markers for with max(`weights`) is
#' maximal, subject to the constraint that no two adjacent markers may
#' be separated by more than `min_distance`.
#'
#' The computation time for the algorithm grows with like the square
#' of the number of markers, like 1 sec for 10k markers
#' but 30 sec for 50k markers. If the number of markers on a chromosome
#' is greater than `max_batch`, the markers are split into batches and
#' the algorithm applied to each batch with min_distance smaller by a
#' factor `min_distance_mult`, and then merged together for one last pass.
#'
#' @seealso [find_dup_markers()], [drop_markers()]
#'
#' @references Broman KW, Weber JL (1999) Method for constructing
#' confidently ordered linkage maps. Genet Epidemiol 16:337--343
#'
#' @examples
#' # read data
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#'
#' # grab genetic map
#' gmap <- grav2$gmap
#'
#' # subset to markers that are >= 1 cM apart
#' gmap_sub <- reduce_markers(gmap, 1)
#'
#' # drop all of the other markers from the cross
#' markers2keep <- unlist(lapply(gmap_sub, names))
#' grav2_sub <- pull_markers(grav2, markers2keep)
#' @export
reduce_markers <-
    function(map, min_distance=1, weights=NULL,
             max_batch=10000, batch_distance_mult=1, cores=1)
{
    if(is.null(map)) stop("map is NULL")
    if(is.cross2(map))
        stop('Input map is a "cross2" object but should be a genetic map')

    if(!is.list(map) && !is.list(weights)) { # single chromosome
        if(!is.null(weights)) weights <- list("1"=weights)
        return(reduce_markers(list("1"=map), min_distance, weights)[[1]])
    }

    if(is.null(weights))
        weights <- lapply(map, function(a) rep(1, length(a)))

    stopifnot(length(map) == length(weights))
    nmar <- vapply(map, length, 1)
    nwts <- vapply(weights, length, 1)
    if(!all(nmar == nwts))
        stop("Different numbers of markers and weights on chr",
             paste(names(map)[nmar != nwts], collapse=" "))

    if(!is_pos_number(min_distance)) stop("min_distance should be a single positive number")

    for(i in seq(along=map)) {
        if(length(map[[i]]) < 2) next
        map[[i]] <- reduce_markers_onechr(map[[i]], min_distance, weights[[i]],
                                          max_batch, batch_distance_mult, cores)
    }

    map
}

reduce_markers_onechr <-
    function(marker_pos, min_distance=1, weights=NULL,
             max_batch=10000, batch_distance_mult=1, cores=1)
{
    num_mar <- length(marker_pos)
    if(is.null(weights)) {
        weights <- rep(1, num_mar)
    }
    stopifnot(num_mar == length(weights))

    cores <- setup_cluster(cores)

    # repeat the following until number hasn't decreased and
    while(num_mar > max_batch) {
        # split markers into groups of length
        batches <- batch_vec(seq_len(num_mar), max_batch, n_cores=1)

#        message(num_mar, " in ", length(batches), " batches")

        batch_func <- function(batch) {
            result <- .reduce_markers(marker_pos[batch], min_distance/batch_distance_mult, weights[batch])
            batch[result]
        }

        sub_map <- unlist( cluster_lapply(cores, batches, batch_func) )

        marker_pos <- marker_pos[sub_map]
        weights <- weights[sub_map]

        if(length(marker_pos) == num_mar) break # there's been no change in number of markers
        num_mar <- length(marker_pos)
    }

    marker_pos[.reduce_markers(marker_pos, min_distance, weights)]
}
