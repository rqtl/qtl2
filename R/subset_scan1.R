# subsetting scan1 output

#' Subset scan1 output
#'
#' Subset the output of [scan1()] by chromosome or column
#'
#' @md
#'
#' @param x An object of class `"scan1"` as returned by
#' [scan1()].
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#' @param chr Vector of chromosomes.
#' @param lodcolumn Vector of integers or character strings indicating the LOD
#' score columns, either as a numeric indexes or column names.
#' @param ... Ignored
#'
#' @export
#' @return Object of same form as input, but subset by chromosome and/or column.
#'
#' @examples
#' \dontrun{
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # pull out chromosome 8
#' out_c8 <- subset(out, map, chr="8")
#'
#' # just the second column on chromosome 2
#' out_c2_spleen <- subset(out, map, "2", "spleen")
#'
#' # all positions, but just the "liver" column
#' out_spleen <- subset(out, map, lodcolumn="spleen")
#' }
subset_scan1 <-
    function(x, map=NULL, chr=NULL, lodcolumn=NULL, ...)
{
    # align scan1 output and map
    if(!is.null(map)) {
        tmp <- align_scan1_map(x, map)
        x <- tmp$scan1
        map <- tmp$map
    }

    x_attr <- attributes(x)
    x_attrnam <- names(x_attr)
    x_class <- class(x)

    # subset to markers within map
    if(!is.null(map)) {
        map_mnames <- map2markernames(map)
        x_mnames <- rownames(x)

        if(!all(x_mnames %in% map_mnames)) {
            x <- x[x_mnames %in% map_mnames,,drop=FALSE]
            x_mnames <- rownames(x)
        }
        if(!all(map_mnames %in% x_mnames)) {
            for(i in seq_along(map)) {
                mn <- names(map[[i]])
                map[[i]] <- map[[i]][mn %in% x_mnames]
            }
        }
    }

    # subset by chromosome
    if(!is.null(chr)) {
        if(is.null(map))
            stop("To subset by chromosome, you need to provide map object.")

        chr <- subset_chr(chr, names(map))

        # rows to keep
        row <- map2chr(map) %in% chr

        x <- x[row,,drop=FALSE]

        # attributes to subset by row
        for(obj in c("SE", "hsq")) {
            if(obj %in% x_attrnam) {
                if(all(chr %in% rownames(x_attr[[obj]])))
                    x_attr[[obj]] <- x_attr[[obj]][chr,,drop=FALSE]
            }
        }

    }

    # subset by column
    if(!is.null(lodcolumn)) {
        if(is.character(lodcolumn)) {
            cols <- colnames(x)

            tmp <- match(lodcolumn, cols)
            if(any(is.na(tmp)))
                stop("Some lodcolumn not found: ", paste(lodcolumn[is.na(tmp)], collapse=", "))
            lodcolumn <- tmp
        }
        if(is.logical(lodcolumn) && length(lodcolumn) != length(cols))
            stop("lodcolumn is logical but not the correct length (", length(cols), ")")

        x <- x[,lodcolumn,drop=FALSE]

        # attributes to subset list
        sublist <- c("sample_size")
        for(obj in sublist) {
            if(obj %in% x_attrnam)
                x_attr[[obj]] <- x_attr[[obj]][lodcolumn]
        }

        # attributes to subset by column
        for(obj in c("SE", "hsq")) {
            if(obj %in% x_attrnam) {
                x_attr[[obj]] <- x_attr[[obj]][,lodcolumn,drop=FALSE]
            }
        }
    }

    # restore attributes
    for(obj in c("SE", "hsq", "sample_size"))
        attr(x, obj) <- x_attr[[obj]]
    class(x) <- x_class

    x
}

#' @export
#' @rdname subset_scan1
subset.scan1 <- subset_scan1


# grab marker names as a vector
map2markernames <-
    function(map)
{
    nam <- unlist(lapply(map, names))
    names(nam) <- NULL
    nam
}

# grab chromosome IDs as a vector
map2chr <-
    function(map)
{
    chr <- rep(names(map), vapply(map, length, 0))
    names(chr) <- map2markernames(map)
    chr
}

# grab positions as a vector
map2pos <-
    function(map)
{
    pos <- unlist(map)
    names(pos) <- map2markernames(map)
    pos
}
