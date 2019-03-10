#' Find gaps in a genetic map
#'
#' Find gaps between markers in a genetic map
#'
#' @param map Genetic map as a list of vectors (each vector is a
#' chromosome and contains the marker positions).
#' @param min_gap Minimum gap length to return.
#'
#' @return
#' Data frame with 6 columns: chromosome, marker to left of gap,
#' numeric index of marker to left, marker to right of gap, numeric
#' index of marker to right, and the length of the gap.
#'
#' @seealso [reduce_map_gaps()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' find_map_gaps(iron$gmap, 40)
#'
#' @export
find_map_gaps <-
    function(map, min_gap=50)
{
    if(is.null(map)) stop("map is NULL")
    if(!is_nonneg_number(min_gap)) stop("min_gap should be a single non-negative number")

    gaps <- NULL
    for(i in seq_along(map)) {
        # interval widths
        d <- diff(map[[i]])

        index <- seq_along(map[[i]])
        lefti <- index[-length(index)]
        righti <- index[-1]
        marnam <- names(map[[i]])

        df <- data.frame(chr=names(map)[i],
                         left_marker=marnam[lefti],
                         left_index=lefti,
                         right_marker=marnam[righti],
                         right_index=righti,
                         gap=d,
                         stringsAsFactors=FALSE)
        df <- df[df$gap >= min_gap, , drop=FALSE]
        gaps <- rbind(gaps, df)
    }

    rownames(gaps) <- seq_len(nrow(gaps))

    gaps
}
