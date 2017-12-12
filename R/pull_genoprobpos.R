#' Pull genotype probabilities for a particular position
#'
#' Pull out the genotype probabilities for a particular position (by
#' name)
#'
#' @md
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param marker A single character string with the name of the
#' position to pull out.
#'
#' @return A matrix of genotype probabilities for the specified position.
#'
#' @seealso [find_marker()], [fit1()]
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
#' @export
pull_genoprobpos <-
    function(genoprobs, marker)
{
    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(length(marker) == 0) stop("marker has length 0")
    if(length(marker) > 1) {
        marker <- marker[1]
        warning("marker should have length 1; using the first value")
    }

    # pull out marker names
    markers <- lapply(genoprobs, function(a) dimnames(a)[[3]])
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
