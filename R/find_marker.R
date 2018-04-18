# find_marker
#' Find markers by chromosome position
#'
#' Find markers closest to specified set of positions, or within a specified interval.
#'
#' @md
#'
#' @param map A map object: a list (corresponding to chromosomes) of
#'     vectors of marker positions. Can also be a snpinfo object (data
#'     frame with columns `chr` and `pos`; marker names taken from
#'     column `snp` or if that doesn't exist from the row names)
#' @param chr A vector of chromosomes
#' @param pos A vector of positions
#' @param interval A pair of positions (provide either `pos` or `interval` but not both)
#'
#' @return A vector of marker names, either closest to the positions
#' specified by `pos`, or within the interval defined by
#' `interval`.
#'
#' @details
#' If `pos` is provided, `interval` should not be, and vice versa.
#'
#' If `pos` is provided, then `chr` and `pos` should
#' either be the same length, or one of them should have length 1 (to
#' be expanded to the length of the other).
#'
#' If `interval` is provided, then `chr` should have length 1.
#'
#' @seealso [find_markerpos()], [find_index_snp()], [pull_genoprobpos()], [pull_genoprobint()]
#' @export
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
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
    function(map, chr, pos=NULL, interval=NULL)
{
    if(is.null(map)) stop("map is NULL")
    if(is.data.frame(map)) { # snp info table; convert to list
        if(!all(c("chr", "pos") %in% colnames(map)))
            stop('If map is a data frame, it should include columns "chr", and "pos".')
        map_chr <- factor(map$chr, unique(map$chr))
        map_pos <- split(map$pos, map_chr)
        mar <- map$snp
        if(is.null(mar)) mar <- rownames(map)
        mar <- split(mar, map_chr)
        for(i in seq_along(map_pos))
            names(map_pos[[i]]) <- mar[[i]]
        map <- map_pos
    }
    chr <- as.character(chr) # treat as a character string

    if(!is.null(interval)) {
        if(!is.null(pos))
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
