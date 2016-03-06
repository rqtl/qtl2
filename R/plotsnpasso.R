#' Plot SNP associations
#'
#' Plot SNP associations, with possible expansion from distinct snps to all snps.
#'
#' @param scan1output Output of \code{\link[qtl2scan]{scan1}} or
#' \code{\link[qtl2scan]{scan1_lmm}}. Should contain an attribute,
#' \code{"snpinfo"}, as when \code{\link[qtl2scan]{scan1}} or
#' \code{\link[qtl2scan]{scan1_lmm}} are run with SNP probabilities
#' produced by \code{\link[qtl2scan]{genoprob_to_snpprob}}.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param cex Character expansion for the points (default 0.5)
#'
#' @param pch Plotting character for the points (default 16)
#'
#' @param drop.hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#'
#' @param col.hilit Color of highlighted points
#'
#' @param col Color of other points
#'
#' @param gap Gap between chromosomes.
#'
#' @param ylim y-axis limits
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @export
#'
plotsnpasso <-
    function(scan1output, show_all_snps=TRUE, drop.hilit=NA,
             col.hilit="violetred", col="black",
             pch=16, cex=0.5, ylim=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    if(show_all_snps)
        scan1output <- expand_snp_results(scan1output)

    if(is.null(ylim))
        ylim <- c(0, max(scan1output, na.rm=TRUE)*1.02)

    if(!is.na(drop.hilit) && !is.null(drop.hilit))
        col <- c(col, col.hilit)[(scan1output >= max(scan1output)-drop.hilit)+1]

    plotscan1(scan1output, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
              gap=gap, add=add, col = col, type="p", cex=cex, pch=pch)
}


# expand snp association results according to snpinfo
expand_snp_results <-
    function(snp_results)
{
    snpinfo <- attr(snp_results, "snpinfo")
    map <- attr(snp_results, "map")

    if(is.null(snpinfo)) stop("No snpinfo found")
    if(is.null(map)) stop("No map found")

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(length(snp_results) != length(unlist(map)))
        stop("length(snp_results) [", length(snp_results), "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    lodindex <- split(seq(along=snp_results), rep(names(map), vapply(map, length, 0)))

    result <- NULL
    for(i in seq(along=map)) {
        map[[i]] <- map[[i]][snpinfo[[i]]$index]
        result <- c(result, snp_results[ lodindex[[i]][ snpinfo[[i]]$index ] ])
    }
    result <- cbind(result)
    attr(result, "map") <- map

    result
}
