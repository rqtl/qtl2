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
#' @seealso \code{\link{plot_scan1}}, \code{\link{plot_coef}}, \code{\link{plot_coefCC}}
#' @export
#'
plot_snpasso <-
    function(scan1output, show_all_snps=TRUE, drop.hilit=NA,
             col.hilit="violetred", col="darkslateblue",
             pch=16, cex=0.5, ylim=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    if(show_all_snps)
        scan1output <- expand_snp_results(scan1output)

    # maximum LOD
    maxlod <- max(scan1output$lod[,1], na.rm=TRUE)

    if(is.null(ylim))
        ylim <- c(0, maxlod*1.02)

    if(!is.na(drop.hilit) && !is.null(drop.hilit))
        col <- c(col, col.hilit)[(scan1output$lod >= maxlod-drop.hilit)+1]

    plot_scan1(scan1output, lodcolumn=1, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
               gap=gap, add=add, col = col, type="p", cex=cex, pch=pch)
}


# expand snp association results according to snpinfo
expand_snp_results <-
    function(snp_results)
{
    snpinfo <- snp_results$snpinfo
    if(is.null(snpinfo)) stop("No snpinfo found")
    map <- snp_results$map
    if(is.null(map)) stop("No map found")

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(length(snp_results$lod) != length(unlist(map)))
        stop("length(snp_results$lod) [", length(snp_results$lod), "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    lodindex <- split(seq(along=snp_results$lod), rep(names(map), vapply(map, length, 0)))

    result <- NULL
    for(i in seq(along=map)) {
        map[[i]] <- snpinfo[[i]]$pos
        names(map[[i]]) <- snpinfo[[i]]$snp
        result <- c(result, snp_results[ lodindex[[i]][ snpinfo[[i]]$index ] ])
    }

    list(lod=result,
         map=map)
}
