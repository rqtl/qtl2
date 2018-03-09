#' Compare two marker maps
#'
#' Compare two marker maps, identifying markers that are only in one
#' of the two maps, or that are in different orders on the two maps.
#'
#' @md
#'
#' @param map1 A list of numeric vectors; each vector gives marker
#' positions for a single chromosome.
#' @param map2 A second map, in the same format as `map1`.
#'
#' @return
#' A data frame containing
#' * `marker` - marker name
#' * `chr_map1` - chromosome ID on `map1`
#' * `pos_map1` - position on `map1`
#' * `chr_map2` - chromosome ID on `map2`
#' * `pos_map2` - position on `map2`
#'
#' @examples
#' # load some data
#' iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
#' gmap <- iron$gmap
#' pmap <- iron$pmap
#'
#' # omit a marker from each map
#' gmap[[7]] <- gmap[[7]][-3]
#' pmap[[8]] <- pmap[[8]][-7]
#' # swap order of a couple of markers on the physical map
#' names(pmap[[9]])[3:4] <- names(pmap[[9]])[4:3]
#' # move a marker to a different chromosome
#' pmap[[10]] <- c(pmap[[10]], pmap[[1]][2])[c(1,3,2)]
#' pmap[[1]] <- pmap[[1]][-2]
#'
#' # compare these messed-up maps
#' compare_maps(gmap, pmap)
#'
#' @importFrom stats setNames
#' @export
compare_maps <-
    function(map1, map2)
{
    if(!is.list(map1) || !is.list(map2)) stop("map1 and map2 should be lists")

    # find markers that are not on both maps
    mar1 <- unlist(lapply(map1, names))
    mar2 <- unlist(lapply(map2, names))
    only1 <- mar1[!(mar1 %in% mar2)]
    only2 <- mar2[!(mar2 %in% mar1)]
    both <- mar1[mar1 %in% mar2]

    # vectors with chr ID and position
    chr1 <- setNames(rep(names(map1), vapply(map1, length, 0)), mar1)
    chr2 <- setNames(rep(names(map2), vapply(map2, length, 0)), mar2)
    pos1 <- setNames(unlist(map1), mar1)
    pos2 <- setNames(unlist(map2), mar2)

    result <- data.frame(marker=character(0),
                         chr_map1=character(0),
                         pos_map1=numeric(0),
                         chr_map2=character(0),
                         pos_map2=numeric(0),
                         stringsAsFactors=FALSE)

    if(length(only1) > 0) {
        result <- rbind(result,
                        data.frame(marker=only1,
                                   chr_map1=chr1[only1],
                                   pos_map1=pos1[only1],
                                   chr_map2=rep(NA, length(only1)),
                                   pos_map2=rep(NA, length(only1)),
                                   stringsAsFactors=FALSE))
    }
    if(length(only2) > 0) {
        result <- rbind(result,
                        data.frame(marker=only2,
                                   chr_map1=rep(NA, length(only2)),
                                   pos_map1=rep(NA, length(only2)),
                                   chr_map2=chr2[only2],
                                   pos_map2=pos2[only2],
                                   stringsAsFactors=FALSE))
    }

    # reduce to common markers
    for(i in seq_along(map1))
        map1[[i]] <- map1[[i]][names(map1[[i]]) %in% both]
    for(i in seq_along(map2))
        map2[[i]] <- map2[[i]][names(map2[[i]]) %in% both]

    # reduce chr and position vectors
    chr1 <- chr1[both]
    chr2 <- chr2[both]
    pos1 <- pos1[both]
    pos2 <- pos2[both]

    # markers on different chromosomes
    diff_chr <- names(chr1)[chr1 != chr2[names(chr1)]]
    if(length(diff_chr) > 0) {
        result <- rbind(result,
                        data.frame(marker=diff_chr,
                                   chr_map1=chr1[diff_chr],
                                   pos_map1=pos1[diff_chr],
                                   chr_map2=chr2[diff_chr],
                                   pos_map2=pos2[diff_chr],
                                   stringsAsFactors=FALSE))
    }

    # reduce to common markers
    same_chr <- both[!(both %in% diff_chr)]
    for(i in seq_along(map1))
        map1[[i]] <- map1[[i]][names(map1[[i]]) %in% same_chr]
    for(i in seq_along(map2))
        map2[[i]] <- map2[[i]][names(map2[[i]]) %in% same_chr]

    # reduce chr and position vectors
    chr1 <- chr1[same_chr]
    chr2 <- chr2[same_chr]
    pos1 <- pos1[same_chr]
    pos2 <- pos2[same_chr]

    # order markers by position on first chromosome
    map2 <- map2[names(map1)]
    for(i in seq_along(map1)) {
        map1[[i]] <- map1[[i]][order(map1[[i]])]
        map2[[i]] <- map2[[i]][names(map1[[i]])]
    }

    # order problems
    for(i in seq_along(map2)) {
        wh <- which(diff(map2[[i]]) < 0)
        wh <- sort(unique(c(wh, wh+1)))
        problems <- names(map1[[i]])[wh]
        if(length(problems) > 0) {
            result <- rbind(result,
                            data.frame(marker=problems,
                                       chr_map1=chr1[problems],
                                       pos_map1=pos1[problems],
                                       chr_map2=chr2[problems],
                                       pos_map2=pos2[problems],
                                       stringsAsFactors=FALSE))
        }
    }

    rownames(result) <- seq_len(nrow(result))

    result
}
