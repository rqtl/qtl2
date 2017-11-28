# find_markerpos
#' Find positions of markers
#'
#' Find positions of markers within a cross object
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param markers A vector of marker names.
#' @param na.rm If TRUE, don't include not-found markers in the
#' results (but issue a warning if some markers weren't found). If
#' FALSE, include those markers with `NA` for chr and position.
#'
#' @return A data frame with chromosome and genetic and physical
#' positions (in columns `"gmap"` and `"pmap"`), with
#' markers as row names.
#'
#' @export
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # find markers
#' find_markerpos(iron, c("D8Mit294", "D11Mit101"))
find_markerpos <-
    function(cross, markers, na.rm=TRUE)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    if(!any(c("gmap", "pmap") %in% names(cross)))
        stop('cross contains neither "gmap" nor "pmap"')

    chr <- rep(NA, length(markers))
    gmap <- rep(NA, length(markers))
    pmap <- rep(NA, length(markers))

    for(i in seq(along=cross$geno)) {
        mn <- colnames(cross$geno[[i]])
        on_this_chr <- (markers %in% mn)
        if(any(on_this_chr)) {
            chr[on_this_chr] <- names(cross$geno)[i]

            if("gmap" %in% names(cross))
                gmap[on_this_chr] <- cross$gmap[[i]][markers[on_this_chr]]

            if("pmap" %in% names(cross))
                pmap[on_this_chr] <- cross$pmap[[i]][markers[on_this_chr]]
        }
    }

    if(na.rm) {
        if(any(is.na(chr)))
            warning("Some markers not found: ", paste(markers[is.na(chr)], collapse=", "))
        gmap <- gmap[!is.na(chr)]
        pmap <- pmap[!is.na(chr)]
        markers <- markers[!is.na(chr)]
        chr <- chr[!is.na(chr)]
    }

    result <- data.frame(chr=chr, gmap=gmap, pmap=pmap,
                         stringsAsFactors=FALSE)
    rownames(result) <- markers
    if(!("pmap" %in% names(cross)))
        result <- result[,-3,drop=FALSE]
    if(!("gmap" %in% names(cross)))
        result <- result[,-2,drop=FALSE]

    result
}
