# convert2cross2
#' Convert R/qtl cross object to new format
#'
#' Convert a cross object from the R/qtl format to the R/qtl2 format
#'
#' @param cross An object of class \code{"cross"}; see
#' \code{\link[qtl]{read.cross}} for details.
#'
#' @return Object of class \code{"cross2"}. For details, see the
#' R/qtl2 Developer Guide. **FIX ME (put the url here)**
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{read_cross2}}
#' @examples
#' library(qtl)
#' data(hyper)
#' hyper2 <- convert2cross2(hyper)
convert2cross2 <-
function(cross)
{
    crosstype <- class(cross)[1]
    check_crosstype(crosstype)
    result <- list(crosstype=crosstype)
    n.ind <- qtl::nind(cross)

    # genetic map, and grab chrtype
    result$gmap <- qtl::pull.map(cross)
    class(result$gmap) <- "list"
    result$is_x_chr <- vapply(result$gmap, function(a) class(a)=="X", TRUE)
    for(i in seq(along=result$gmap))
        class(result$gmap[[i]]) <- "numeric"

    # ids
    ids <- qtl::getid(cross)
    if(is.null(ids)) ids <- as.character(1:n.ind)

    # split out genotype data
    result$geno <- lapply(cross$geno, "[[", "data")
    for(i in seq(along=result$geno))
        rownames(result$geno[[i]])  <- ids

    # sex/pgm
    sexpgm <- qtl::getsex(cross)
    if(is.null(sexpgm$sex))
        result$sex <- rep(FALSE, n.ind)
    else result$sex <- sexpgm$sex
    names(result$sex) <- ids

    if(is.null(sexpgm$pgm))
        result$cross_info <- matrix(ncol=0, nrow=n.ind)
    else result$cross_info <- matrix(sexpgm$pgm)
    rownames(result$cross_info) <- ids

    # convert X chr genotypes
    if(any(result$is_x_chr) && crosstype %in% c("bc", "f2", "bcsft")) { # bcsft not really supported yet
        for(i in which(result$is_x_chr)) {
            result$geno[[i]] <- qtl::reviseXdata(crosstype, "simple",
                                                 sexpgm, geno=result$geno[[i]],
                                                 cross.attr=attributes(cross), force=TRUE)
        }
    }

    # in genotypes, replace NAs with 0s
    for(i in seq(along=result$geno))
        result$geno[[i]][is.na(result$geno[[i]])] <- 0

    # phenotypes: pull out numeric columns and put the rest in covariates
    phe <- cross$pheno

    # put sex and pgm in covar data even if numeric
    phe.names <- colnames(phe)
    is_sex <- is_pgm <- rep(FALSE, ncol(phe))
    is_sex[grep("^[Ss][Ee][Xx]$", phe.names)] <- TRUE
    is_pgm[grep("^[Pp][Gg][Mm]$", phe.names)] <- TRUE
    numercol <- vapply(phe, is.numeric, TRUE) & (!is_sex) & (!is_pgm)
    covar <- phe[,!numercol,drop=FALSE]
    if(ncol(covar) > 0) {
        result$covar <- covar
        rownames(result$covar) <- ids
    }
    phe <- phe[, numercol, drop=FALSE]
    if(ncol(phe) > 0) {
        result$pheno <- as.matrix(phe)
        rownames(result$pheno) <- ids
        storage.mode(result$pheno) <- "double"
    }

    class(result) <- "cross2"
    result
}
