# find_marker
#' Find markers by chromosome position
#'
#' Find markers closest to specified set of positions, or within a specified interval.
#'
#' @param map A map object: a list (corresponding to chromosomes) of vectors of marker positions.
#' @param chr A vector of chromosomes
#' @param pos A vector of positions
#' @param interval A pair of positions (provide either \code{pos} or \code{interval} but not both)
#'
#' @return A vector of marker names, either closest to the positions
#' specified by \code{pos}, or within the interval defined by
#' \code{interval}.
#'
#' @details
#' If \code{pos} is provided, \code{interval} should not be, and vice versa.
#'
#' If \code{pos} is provided, then \code{chr} and \code{pos} should
#' either be the same length, or one of them should have length 1 (to
#' be expanded to the length of the other).
#'
#' If \code{interval} is provided, then \code{chr} should have length 1.
#'
#' @export
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # find markers by their genetic map positions
#' find_marker(iron$gmap, c(8, 11), c(37.7, 56.9))
#'
#' # find markers by their physical map positions (two markers on chr 7)
#' find_marker(iron$pmap, 7, c(44.2, 108.9))
#'
#' # find markers in an interval
#' find_marker(iron$pmap, 16, interval=c(35, 80))
find_marker <-
    function(map, chr, pos, interval=NULL)
{
    chr <- as.character(chr) # treat as a character string

    if(!is.null(interval)) {
        if(!missing(pos) && !is.null(pos))
            stop("Only one of pos and interval should be given")
        if(length(chr) > 1)
            stop("If interval is provided, chr should have length 1")
        if(length(interval) != 2)
            stop("interval should have length 2")

        # make sure it's sorted
        interval <- sort(interval)

        if(!(chr %in% names(map)))
            stop("Can't find chromosome ", chr)

        map <- map[[chr]]
        return(names(map)[map >= interval[1] & map <= interval[2]])
    }
    else {
        if(length(chr) != length(pos)) {
            if(length(chr) > 1 && length(pos) > 1)
                stop("chr and pos should be the same length, or at least one of length 1")
            if(length(chr)==1 && length(pos) > 1)
                chr <- rep(chr, length(pos))
            if(length(pos)==1 && length(chr) > 1)
                pos <- rep(pos, length(chr))
        }
        if(length(chr) > 1) {
            found <- chr %in% names(map)
            if(!all(found)) {
                if(!any(found))
                    stop("None of the chromosomes found in the map")
                warning("Some chr not found: ", paste(chr[!found], collapse=", "))
                chr <- chr[found]
                pos <- pos[found]
            }

            result <- rep("", length(chr))
            for(i in seq(along=chr))
                result[i] <- find_marker(map, chr[i], pos[i])
            return(result)
        }

        # one position
        if(!(chr %in% names(map)))
            stop("chr ", chr, " not found")
        map <- map[[chr]]
        d <- abs(map - pos)
        closest <- which(d==min(d))
        if(length(closest) > 1) closest <- sample(closest, 1)
        return(names(map)[closest])
    }
}
