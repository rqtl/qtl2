#' Interpolate between maps
#'
#' Use interpolate to convert from one map to another
#'
#' @param map The map to be interpolated; a list of vectors.
#' @param oldmap Map with positions in the original scale, as in \code{map_to_change}.
#' @param newmap Map with positions in the new scale.
#'
#' @export
#' @return Object of same form as input \code{map} but in the units as in \code{newmap}.
interp_map <-
    function(map, oldmap, newmap)
{
    chr <- names(map)
    result <- vector("list", length(map))
    names(result) <- chr

    chr_found <- (chr %in% names(oldmap)) & (chr %in% names(newmap))
    if(!all(chr_found)) {
        stop("Chromosomes not found in both old and new maps: ",
             paste(chr[!chr_found], collapse=", "))
    }

    for(thechr in chr) {
        om <- oldmap[[thechr]]
        nm <- newmap[[thechr]]
        if(length(om) != length(nm) || !all(names(om) == names(nm)))
            stop("Old and new maps differ on chr ", thechr)
        result[[thechr]] <- interpolate_map(map[[thechr]], om, nm)
    }

    result
}
