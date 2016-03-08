# convert2scanone
#
#' Convert scan1 results to the scanone format
#'
#' Convert the results of \code{\link{scan1}} or
#' \code{\link{scan1_lmm}} to the form used by the R/qtl function
#' \code{scanone}.
#'
#' @param output Matrix of LOD scores, as calculated by
#' \code{\link{scan1}} or \code{\link{scan1_lmm}}.
#'
#' @return A data frame with class \code{"scanone"}, containing
#' chromosome and position columns followed by the LOD scores in
#' \code{output}.
#'
#' @examples
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#' out <- scan1(probs, pheno, covar, Xcovar)
#'
#' out_rev <- convert2scanone(out)
#'
#' @export
convert2scanone <-
    function(output)
{
    map <- output$map
    n <- sapply(map, length)

    # to handle output of either scan1() or scan1coef()
    if("lod" %in% names(output))
        lod <- output$lod
    else if("coef" %in% names(output))
        lod <- output$coef
    else stop("Neither lod nor coef found.")

    stopifnot(sum(n) == nrow(lod))

    out <- data.frame(chr=factor(rep(names(map), n), levels=names(map)),
                      pos=unlist(map),
                      as.data.frame(lod))
    rownames(out) <- rownames(output)

    class(out) <- c("scanone", "data.frame")
    out
}
