#' Plot SNP associations
#'
#' Plot SNP associations, with possible expansion from distinct snps to all snps.
#'
#' @param scan1output Output of \code{\link[qtl2scan]{scan1}}.  Should
#' contain an attribute, \code{"snpinfo"}, as when
#' \code{\link[qtl2scan]{scan1}} are run with SNP probabilities
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
#' @examples
#' \dontrun{
#' # load example DO data from web
#' library(qtl2geno)
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#'
#' # subset to chr 2
#' DOex <- DOex[,"2"]
#'
#' # calculate genotype probabilities and convert to allele probabilities
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' apr <- genoprob_to_alleleprob(pr)
#'
#' # download snp info from web
#' tmpfile <- tempfile()
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/c2_snpinfo.rds")
#' download.file(file, tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#'
#' # calculate strain distribution patterns
#' snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
#'
#' # switch map in allele probabilities to Mbp
#' apr$map <- DOex$pmap
#'
#' # convert to snp probabilities
#' snppr <- genoprob_to_snpprob(apr, snpinfo)
#'
#' # perform SNP association analysis (here, ignoring residual kinship)
#' library(qtl2scan)
#' out_snps <- scan1(snppr, DOex$pheno)
#'
#' # plot results
#' library(qtl2plot)
#' plot_snpasso(out_snps)
#'
#' # can also just type plot()
#' plot(out_snps)
#'
#' # plot just subset of distinct SNPs
#' plot_snpasso(out_snps, show_all_snps=FALSE)
#'
#' # highlight the top snps (with LOD within 1.5 of max)
#' plot(out_snps, drop.hilit=1.5)
#' }
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
               gap=gap, add=add, col = col, type="p", cex=cex, pch=pch, ...)
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
        result <- rbind(result, snp_results$lod[lodindex[[i]],,drop=FALSE][snpinfo[[i]]$index,,drop=FALSE])
    }

    list(lod=result,
         map=map)
}
