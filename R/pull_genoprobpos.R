#' Pull genotype probabilities for a particular position
#'
#' Pull out the genotype probabilities for a particular position (by
#' name)
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param map A map object: a list (corresponding to chromosomes) of
#'     vectors of marker positions. Can also be a snpinfo object (data
#'     frame with columns `chr` and `pos`; marker names taken from
#'     column `snp` or if that doesn't exist from the row names)
#' @param chr A chromosome ID
#' @param pos A numeric position
#' @param marker A single character string with the name of the
#' position to pull out.
#'
#' @return A matrix of genotype probabilities for the specified position.
#'
#' @details Provide either a marker/pseudomarker name (with the argument `marker`)
#' or all of `map`, `chr`, and `pos`.
#'
#' @seealso [find_marker()], [fit1()], [pull_genoprobint()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,"8"]}
#' gmap <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, gmap, error_prob=0.002)
#'
#' pmar <- find_marker(gmap, 8, 40)
#' pr_8_40 <- pull_genoprobpos(pr, pmar)
#'
#' pr_8_40_alt <- pull_genoprobpos(pr, gmap, 8, 40)
#'
#' @export
pull_genoprobpos <-
    function(genoprobs, map=NULL, chr=NULL, pos=NULL, marker=NULL)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")

    if(is.character(map) && is.null(chr) && is.null(pos) && is.null(marker)) {
        # treat this as map -> marker
        marker <- map
        map <- NULL
    }

    if(!is.null(marker) && (!is.null(map) || !is.null(chr) || !is.null(pos))) {
        warning("Provide either {map,chr,pos} or marker, not both; using marker")
    }
    if(is.null(marker)) {
        marker <- find_marker(map, chr, pos)
    }

    if(length(marker) == 0) stop("marker has length 0")
    if(length(marker) > 1) {
        marker <- marker[1]
        warning("marker should have length 1; using the first value")
    }

    # pull out marker names
    markers <- dimnames(genoprobs)[[3]]
    # vector of chromosomes
    chr <- rep(names(genoprobs), lapply(markers, length))
    # vector of indexes
    index <- unlist(lapply(markers, function(a) seq_along(a)))
    # markers into single vector
    markers <- unlist(markers)

    # identify the particular position
    wh <- (marker == markers)

    # check that it's there exactly onces
    n_found <- sum(wh)
    if(n_found==0) stop('marker "', marker, '" not found')
    if(n_found > 1) stop('marker "', marker, '" appears ', n_found, ' times')

    # pull out that set of probabilities
    result <- genoprobs[[chr[wh]]][,,index[wh], drop=FALSE]

    # make sure it's a plain matrix
    d <- dim(result)
    dn <- dimnames(result)
    result <- matrix(result[,,1], nrow=d[1], ncol=d[2])
    dimnames(result) <- dn[1:2]

    result
}
