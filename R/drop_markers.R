# drop_markers
#' Drop markers from a cross2 object
#'
#' Drop a vector of markers from a cross2 object.
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param markers A vector of marker names.
#'
#' @return The input `cross` with the specified markers removed.
#'
#' @export
#' @seealso [pull_markers()], [drop_nullmarkers()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' markers2drop <- c("BH.342C/347L-Col", "GH.94L", "EG.357C/359L-Col", "CD.245L", "ANL2")
#' grav2_rev <- drop_markers(grav2, markers2drop)
drop_markers <-
    function(cross, markers)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    # to keep track of which markers were found
    found <- rep(FALSE, length(markers))
    names(found) <- markers

    chr2drop <- NULL
    for(i in seq(along=cross$geno)) {
        mn <- colnames(cross$geno[[i]])
        to_drop <- mn %in% markers

        if(any(to_drop)) {
            cross$geno[[i]] <- cross$geno[[i]][,!to_drop,drop=FALSE]
            found[mn[to_drop]] <- TRUE

            if("founder_geno" %in% names(cross))
                cross$founder_geno[[i]] <- cross$founder_geno[[i]][,!to_drop,drop=FALSE]
            if("gmap" %in% names(cross))
                cross$gmap[[i]] <- cross$gmap[[i]][!to_drop]
            if("pmap" %in% names(cross))
                cross$pmap[[i]] <- cross$pmap[[i]][!to_drop]

            # nothing left?
            if(all(to_drop)) chr2drop <- c(chr2drop, i)
        }
    }

    if(length(chr2drop) > 0) {
        cross$geno <- cross$geno[-chr2drop]

        if("founder_geno" %in% names(cross))
            cross$founder_geno <- cross$founder_geno[-chr2drop]
        if("gmap" %in% names(cross))
            cross$gmap <- cross$gmap[-chr2drop]
        if("pmap" %in% names(cross))
            cross$pmap <- cross$pmap[-chr2drop]
        if("is_x_chr" %in% names(cross))
            cross$is_x_chr <- cross$is_x_chr[-chr2drop]
    }
    if(length(cross$geno) == 0)
        stop("No markers left")

    if(any(!found))
        warning("Some markers not found: ", paste(markers[!found], collapse=", "))

    cross
}


# drop_nullmarkers
#' Drop markers with no genotype data
#'
#' Drop markers with no genotype data (or no informative genotypes)
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param quiet If FALSE, print information about how many markers were dropped.
#'
#' @return The input `cross` with the uninformative markers removed.
#'
#' @details We omit any markers that have completely missing data, or
#' if founder genotypes are present (e.g., for Diversity Outbreds),
#' the founder genotypes are missing or are all the same.
#'
#' @export
#' @seealso [drop_markers()], [pull_markers()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' # make a couple of markers missing
#' grav2$geno[[2]][,c(3,25)] <- 0
#' grav2_rev <- drop_nullmarkers(grav2)
drop_nullmarkers <-
    function(cross, quiet=FALSE)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    ngeno <- n_typed(cross, "marker", "count")
    if(any(ngeno==0)) {
        if(!quiet) message("Dropping ", sum(ngeno==0), " markers with no data")
        cross <- drop_markers(cross, names(ngeno)[ngeno==0])
    }

    # look at founder_geno
    if("founder_geno" %in% names(cross)) {
        to_drop <- NULL
        for(i in seq(along=cross$founder_geno)) {
            fg <- cross$founder_geno[[i]]

            # any completely missing?
            ntyped <- colSums(fg!=0)
            if(any(ntyped==0))
                to_drop <- c(to_drop, names(ntyped)[ntyped==0])

            # any with none missing but all the same
            noninf <- apply(fg[,ntyped==nrow(fg),drop=FALSE], 2,
                            function(a) length(unique(a))==1 & all(a != 0))
            if(any(noninf))
                to_drop <- c(to_drop, names(noninf)[noninf])
        }

        if(length(to_drop) > 0) {
            if(!quiet) message("Dropping ", length(to_drop), " noninformative markers")
            cross <- drop_markers(cross, to_drop)
        }
    }

    cross
}


# pull_markers
#' Drop all but a specified set of markers
#'
#' Drop all markers from a cross2 object expect those in a specified vector.
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param markers A vector of marker names.
#'
#' @return The input `cross` with only the specified markers.
#'
#' @export
#' @seealso [drop_markers()], [drop_nullmarkers()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' markers2drop <- c("BH.342C/347L-Col", "GH.94L", "EG.357C/359L-Col", "CD.245L", "ANL2")
#' grav2_rev <- pull_markers(grav2, markers2drop)
pull_markers <-
    function(cross, markers)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    mn <- marker_names(cross)
    to_drop <- !(mn %in% markers)
    not_found <- !(markers %in% mn)

    if(all(to_drop))
        stop("No markers found")

    if(any(not_found))
        warning("Some markers not found: ", paste(markers[not_found], collapse=", "))

    drop_markers(cross, mn[to_drop])
}
