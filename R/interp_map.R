#' Interpolate between maps
#'
#' Use interpolate to convert from one map to another
#'
#' @md
#'
#' @param map The map to be interpolated; a list of vectors.
#' @param oldmap Map with positions in the original scale, as in `map`.
#' @param newmap Map with positions in the new scale.
#'
#' @return Object of same form as input `map` but in the units as in `newmap`.
#'
#' @examples
#' # load example data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
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
    if(is.null(map)) stop("map is NULL")
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

    dropped_old <- dropped_new <- NULL
    for(thechr in chr) {
        om <- oldmap[[thechr]]
        nm <- newmap[[thechr]]
        nom <- names(om)
        nnm <- names(nm)

        tm <- map[[thechr]]
        ntm <- names(tm)

        # markers in oldmap that aren't in newmap?
        if((is.null(nom) || is.null(nnm)) && length(om) != length(nm))
            stop("If length(oldmap) != length(newmap), they both need marker names")

        if(!is.null(nom) && !is.null(nnm)) {
            if(length(unique(nom)) != length(nom) ||
               length(unique(nnm)) != length(nnm))
                stop("marker names should be distinct")

            drop_om <- !(nom %in% nnm)
            if(any(drop_om)) {
                dropped_old <- c(dropped_old, nom[drop_om])
                om <- om[!drop_om]
                nom <- names(om)
            }

            # markers in newmap that aren't in oldmap?
            drop_nm <- !(nnm %in% nom)
            if(any(drop_nm)) {
                dropped_new <- c(dropped_new, nnm[drop_nm])
                nm <- nm[!drop_nm]
                nnm <- names(nm)
            }
        }

        # need at least two markers left, and need them to be in the same order
        if(length(nm) < 2 || length(om) < 2)
            stop("Less than 2 matching markers on chr ", thechr)
        if(length(om) != length(nm) || (!is.null(nom) && !is.null(nnm) && !all(nom == nnm)))
            stop("Old and new maps not aligned on chr ", thechr)

        result[[thechr]] <- interpolate_map(tm, om, nm)
        names(result[[thechr]]) <- ntm

        # make sure markers that are in both map and newmap are kept in the newmap position
        if(!is.null(ntm) && !is.null(nom) && !is.null(nnm) && any(ntm %in% nom)
           && length(unique(ntm)) == length(ntm)) { # skip if marker names aren't distinct
            overlap <- ntm[ntm %in% nom] # markers in both old and to-be-transformed map
            overlap <- overlap[abs(tm[overlap] - om[overlap]) < 1e-6] # ...that are in the same position

            # insert the position in newmap
            result[[thechr]][overlap] <- nm[overlap]
        }
    }

    if(length(dropped_old) > 1)
        warning("Dropped ", length(dropped_old), " markers not present in new map")
    if(length(dropped_new) > 1)
        warning("Dropped ", length(dropped_new), " markers not present in old map")

    result
}
