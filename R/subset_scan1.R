# subsetting scan1 output

#' Subset scan1 output
#'
#' Subset the output of \code{\link{scan1}} by chromosome or column
#'
#' @param x An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}} or \code{\link{scan1_lmm}}.
#' @param chr Vector of chromosomes.
#' @param columns Vector of integers or character strings indicating the LOD
#' score columns, either as a numeric indexes or column names.
#' @param ... Ignored
#'
#' @export
#' @return Object of same form as input, but subset by chromosome and/or column.
#'
#' @examples
#' \dontrun{
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, covar, Xcovar)
#'
#' # pull out chromosome 8
#' out-c8 <- subset(out, chr="8")
#'
#' # just the second column on chromosome 2
#' out_c2_spleen <- subset(out, "2", "spleen")
#' }
subset.scan1 <-
    function(x, chr, columns, ...)
{
    map <- attr(x, "map")
    snpinfo <- attr(x, "snpinfo")
    se <- attr(x, "SE")
    hsq <- attr(x, "hsq")
    sample_size <- attr(x, "sample_size")
    all_attr <- attributes(x)
    lod <- unclass(x)

    # subset by chromosome
    if(!missing(chr) && !is.null(chr)) {
        # selected chromosomes
        chr_found <- chr %in% names(map)
        if(!all(chr_found))
            stop("Not all chr found: ", paste(chr[!chr_found], collapse=", "))

        # rows to keep
        row <- chr_scan1(x) %in% chr

        # subset by chr
        lod <- lod[row,,drop=FALSE]
        map <- map[chr]
        if(!is.null(snpinfo)) snpinfo <- snpinfo[chr]
        if(!is.null(se)) se <- se[row,,drop=FALSE]
        if(!is.null(hsq) && all(chr %in% rownames(hsq)))
            hsq <- hsq[chr,,drop=FALSE]
    }

    if(!missing(columns) && !is.null(columns)) {
        if(is.character(columns)) {
            tmp <- match(columns, colnames(lod))
            if(any(is.na(tmp)))
                stop("Some columns not found: ", paste(columns[is.na(tmp)], collapse=", "))
            columns <- tmp
        }
        if(is.logical(columns) && length(columns) != ncol(lod))
            stop("columns is logical but not the correct length (", ncol(lod), ")")
        lod <- lod[,columns,drop=FALSE]
        if(!is.null(se)) se <- se[,columns,drop=FALSE]
        if(!is.null(hsq)) hsq <- hsq[,columns,drop=FALSE]
        if(!is.null(sample_size)) sample_size <- sample_size[columns]
    }

    # add the attributes back
    attr(lod, "map") <- map
    attr(lod, "snpinfo") <- snpinfo
    attr(lod, "SE") <- se
    attr(lod, "hsq") <- hsq
    attr(lod, "sample_size") <- sample_size

    # all other attributes
    attr2skip <- c("dim", "dimnames", "map", "snpinfo", "hsq", "SE", "sample_size")
    attrnam <- names(all_attr)
    for(i in seq(along=all_attr)) {
        if(attrnam[i] %in% attr2skip) next
        attr(lod, attrnam[i]) <- all_attr[[i]]
    }

    lod
}


# chr_scan1: grab chromosome IDs as a vector
chr_scan1 <-
    function(scan1_output)
{
    map <- attr(scan1_output, "map")

    if(is.null(map))
        stop("No map attribute found.")

    chr <- rep(names(map), vapply(map, length, 0))
    names(chr) <- rownames(scan1_output)
    chr
}


# pos_scan1: grab positions as a vector
pos_scan1 <-
    function(scan1_output)
{
    map <- attr(scan1_output, "map")

    if(is.null(map))
        stop("No map attribute found.")

    pos <- unlist(map)
    names(pos) <- rownames(scan1_output)
    pos
}
