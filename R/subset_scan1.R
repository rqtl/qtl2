# subsetting scan1 output

#' Subset scan1 output
#'
#' Subset the output of \code{\link{scan1}} by chromosome or column
#'
#' @param x An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}}.
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
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # pull out chromosome 8
#' out-c8 <- subset(out, chr="8")
#'
#' # just the second column on chromosome 2
#' out_c2_spleen <- subset(out, "2", "spleen")
#' }
subset_scan1 <-
    function(x, chr=NULL, lodcolumn=NULL, ...)
{
    map <- x$map

    # subset by chromosome
    if(!is.null(chr)) {
        chr <- as.character(chr)

        # selected chromosomes
        chr_found <- chr %in% names(map)
        if(!all(chr_found))
            stop("Not all chr found: ", paste(chr[!chr_found], collapse=", "))

        # rows to keep
        row <- chr_scan1(x) %in% chr

        # objects to subset rows
        subrows <- c("lod", "coef", "SE")
        for(obj in subrows) {
            if(obj %in% names(x))
                x[[obj]] <- x[[obj]][row,,drop=FALSE]
        }

        # objects to subset list
        sublist <- c("map", "snpinfo", "is_x_chr")
        for(obj in sublist) {
            if(obj %in% names(x))
                x[[obj]] <- x[[obj]][chr]
        }

        if(!is.null(x$hsq)) {
            if(all(chr %in% rownames(x$hsq)))
                x$hsq <- x$hsq[chr,,drop=FALSE]
        }
    }

    if(!is.null(lodcolumn)) {
        if(is.character(lodcolumn)) {
            if(!is.null(x$lod)) cols <- colnames(x$lod)
            else if(!is.null(x$coef)) cols <- colnames(x$coef)
            else stop("Neither lod nor coef found.")

            tmp <- match(lodcolumn, cols)
            if(any(is.na(tmp)))
                stop("Some lodcolumn not found: ", paste(lodcolumn[is.na(tmp)], collapse=", "))
            lodcolumn <- tmp
        }
        if(is.logical(lodcolumn) && length(lodcolumn) != length(cols))
            stop("lodcolumn is logical but not the correct length (", length(cols), ")")

        subcol <- c("coef", "lod", "SE", "hsq")
        for(obj in subcol) {
            if(obj %in% names(x))
                x[[obj]] <- x[[obj]][,lodcolumn,drop=FALSE]
        }
        x$sample_size <- x$sample_size[lodcolumn]
    }

    x
}

#' @export
#' @rdname subset_scan1
subset.scan1 <- subset_scan1

#' @export
#' @rdname subset_scan1
`[.scan1` <-
    function(x, chr=NULL, lodcolumn=NULL)
    subset(x, chr, lodcolumn)


# chr_scan1: grab chromosome IDs as a vector
chr_scan1 <-
    function(scan1_output)
{
    map <- scan1_output$map

    if(is.null(map))
        stop("No map found.")

    chr <- rep(names(map), vapply(map, length, 0))

    if("lod" %in% names(scan1_output))
        names(chr) <- rownames(scan1_output$lod)
    else if("coef" %in% names(scan1_output))
        names(chr) <- rownames(scan1_output$coef)
    else stop("Neither lod nor coef found in scan1_output.")

    chr
}


# pos_scan1: grab positions as a vector
pos_scan1 <-
    function(scan1_output)
{
    map <- scan1_output$map

    if(is.null(map))
        stop("No map found.")

    pos <- unlist(map)

    if("lod" %in% names(scan1_output))
        names(pos) <- rownames(scan1_output$lod)
    else if("coef" %in% names(scan1_output))
        names(pos) <- rownames(scan1_output$coef)
    else stop("Neither lod nor coef found in scan1_output.")

    pos
}
