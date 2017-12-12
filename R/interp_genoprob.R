# interp_genoprob
#' Interpolate genotype probabilities
#'
#' Linear interpolation of genotype probabilities, mostly to get two sets onto the same map for comparison purposes.
#'
#' @md
#'
#' @param probs Genotype probabilities, as calculated from
#' [calc_genoprob()].
#' @param map List of vectors of map positions.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return An object like the input `probs` but with additional
#' positions present in `map`.
#'
#' @details We reduce `probs` to the positions present in `map` and then
#' interpolate the genotype probabilities at additional positions
#' in `map` by linear interpolation using the two adjacent
#' positions. Off the ends, we just copy over the first or last
#' value unchanged.
#'
#' In general, it's better to use [insert_pseudomarkers()] and
#' [calc_genoprob()] to get genotype probabilities at additional
#' positions along a chromosome. This function is a **very** crude
#' alternative that was implemented in order to compare genotype
#' probabilities derived by different methods, where we first need to
#' get them onto a common set of positions.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[1:20,c("1", "X")]}
#' probs <- calc_genoprob(iron, iron$gmap, error_prob=0.002)
#'
#' # you generally wouldn't want to do this, but this is an illustration
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' probs_map <- interp_genoprob(probs, map)
#'
#' @seealso [calc_genoprob()]
#'
#' @export
interp_genoprob <-
    function(probs, map, cores=1)
{
    if(is.null(probs)) stop("probs is NULL")
    if(is.null(map)) stop("map is NULL")

    gchr <- names(probs)
    mchr <- names(map)
    if(!all(gchr %in% mchr) || !all(mchr %in% gchr)) {
        cchr <- gchr[gchr %in% mchr]
        probs <- probs[,cchr]
        map <- map[cchr]
    }

    # set up cluster
    cores <- setup_cluster(cores)
    result <- cluster_lapply(cores, seq_along(map), interp_genoprob_onechr, probs=probs, map=map)

    names(result) <- names(probs)
    for(x in c("crosstype", "is_x_chr", "alleles", "alleleprobs", "class"))
        attr(result, x) <- attr(probs, x)

    result
}


interp_genoprob_onechr <-
    function(chr, probs, map)
{
    probs <- probs[[chr]]
    map <- map[[chr]]

    d <- dim(probs)
    markers <- dimnames(probs)[[3]]
    pmar <- names(map)

    keep <- (markers %in% pmar)
    if(!any(keep)) stop("No overlapping markers on chr ", chr)
    if(!all(keep)) {
        probs <- probs[,,keep]
        markers <- markers[keep]
    }

    # check that the marker order didn't change
    m <- match(markers, pmar)
    if(any(diff(m) < 0)) stop("probs positions out of order on chr ", chr)

    pos_index <- match(pmar, markers) - 1 # indexes start at 0
    if(!any(is.na(pos_index))) return(probs) # no new positions
    pos_index[is.na(pos_index)] <- -1

    result <- .interp_genoprob_onechr(probs, map, pos_index)
    dimnames(result) <- c(dimnames(probs)[1:2], list(pmar))

    result
}
