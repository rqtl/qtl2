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
#' @return An object like the input `probs` but with additional positions present in `map`.
#'
#' @details We reduce `probs` to those present in `map` and then
#' interpolate the genotype probabilities at additional positions
#' in `map` by linear interpolation using the two adjacent
#' positions. Off the ends, we just copy over the first or last
#' value unchanged.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
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
    gchr <- names(probs)
    mchr <- names(map)
    if(!all(gchr %in% mchr) || !all(mchr %in% gchr)) {
        cchr <- gchr[gchr %in% mchr]
        probs <- probs[,cchr]
        map <- map[cchr]
    }

    # set up cluster
    cores <- setup_cluster(cores, quiet)
    result <- cluster_lapply(cores, seq_along(map), interp_genoprob_onechr, probs=probs, map=map)

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

    is_new_pos <- !(pmar %in% markers)
    if(!any(is_new_pos)) return(probs)

    result <- .interp_genoprob_onechr(probs, map, is_new_pos)
    dimnames(result) <- c(dimnames(probs)[1:2], list(pmar))

    result
}
