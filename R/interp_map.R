#' Interpolate between maps
#'
#' Use interpolate to convert from one map to another
#'
#' @param map The map to be interpolated; a list of vectors.
#' @param oldmap Map with positions in the original scale, as in \code{map}.
#' @param newmap Map with positions in the new scale.
#'
#' @return Object of same form as input \code{map} but in the units as in \code{newmap}.
#'
#' @examples
#' # load example data
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # positions to interpolate from cM to Mbp
#' tointerp <- list("7" = c(pos7.1= 5, pos7.2=15, pos7.3=25),
#'                  "9" = c(pos9.1=20, pos9.2=40))
#'
#' interp_map(tointerp, iron$gmap, iron$pmap)
#'
#' @export
interp_map <-
    function(map, oldmap, newmap)
{
    if(is.null(oldmap)) stop("oldmap is NULL")
    if(is.null(newmap)) stop("newmap is NULL")

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
