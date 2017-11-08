#' Plot genotype probabilities for one individual on one chromosome.
#'
#' Plot the genotype probabilities for one individual on one chromosome, as a heat map.
#'
#' @param probs Genotype probabilities (as produced by \code{\link[qtl2geno]{calc_genoprob}})
#' or allele dosages (as produced by \code{\link[qtl2geno]{genoprob_to_alleleprob}}).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosome to plot; a single character string.
#' @param geno Optional vector of genotypes or alleles to be shown
#' (vector of integers or character strings)
#' @param color_scheme Color scheme for the heatmap (ignored if \code{col} is provided).
#' @param col Optional vector of colors for the heatmap.
#' @param threshold Threshold for genotype probabilities; only genotypes that achieve
#' this value somewhere on the chromosome will be shown.
#' @param swap_axes If TRUE, swap the axes, so that the genotypes are
#' on the x-axis and the chromosome position is on the y-axis.
#' @param hlines Position of horizontal grid lines (use \code{NA} to avoid lines).
#' @param hlines_col Color of horizontal grid lines.
#' @param hlines_lwd Line width of horizontal grid lines.
#' @param hlines_lty Line type of horizontal grid lines.
#' @param vlines Position of vertical grid lines (use \code{NA} to avoid lines).
#' @param vlines_col Color of vertical grid lines.
#' @param vlines_lwd Line width of vertical grid lines.
#' @param vlines_lty Line type of vertical grid lines.
#' @param ... Additional graphics parameters passed to \code{\link[graphics]{image}}.
#'
#' @examples
#' # load data and calculate genotype probabilities
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' iron <- iron[,"2"] # subset to chr 2
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # plot the probabilities for the individual labeled "262"
#' #  (white = 0, black = 1)
#' plot_genoprob(pr, map, ind="262")
#'
#' # change the x-axis label
#' plot_genoprob(pr, map, ind="262", xlab="Position (cM)")
#'
#' # swap the axes so that the chromosome runs vertically
#' plot_genoprob(pr, map, ind="262", swap_axes=TRUE, ylab="Position (cM)")
#'
#' # This is more interesting for a Diversity Outbred mouse example
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' # subset to chr 2 and X and individuals labeled "232" and "256"
#' DOex <- DOex[c("232", "256"), c("2", "X")]
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' # plot individual "256" on chr 2 (default is to pick first chr in the probs)
#' plot_genoprob(pr, DOex$pmap, ind="256")
#'
#' # omit states that never have probability >= 0.5
#' plot_genoprob(pr, DOex$pmap, ind="256", threshold=0.05)
#'
#' # X chr male 232: just show the AY-HY genotype probabilities
#' plot_genoprob(pr, DOex$pmap, ind="232", chr="X", geno=paste0(LETTERS[1:8], "Y"))
#' # could also indicate genotypes by number
#' plot_genoprob(pr, DOex$pmap, ind="232", chr="X", geno=37:44)
#' # and can use negative indexes
#' plot_genoprob(pr, DOex$pmap, ind="232", chr="X", geno=-(1:36))
#'
#' # X chr female 256: just show the first 36 genotype probabilities
#' plot_genoprob(pr, DOex$pmap, ind="256", chr="X", geno=1:36)
#'
#' # again, can give threshold to omit genotypes whose probabilities never reach that threshold
#' plot_genoprob(pr, DOex$pmap, ind="256", chr="X", geno=1:36, threshold=0.5)
#'
#' # can also look at the allele dosages
#' apr <- genoprob_to_alleleprob(pr)
#' plot_genoprob(apr, DOex$pmap, ind="232")
#' }
#'
#' @export
#' @importFrom graphics image par axis title box
#' @importFrom grDevices gray
plot_genoprob <-
    function(probs, map, ind=1, chr=NULL, geno=NULL,
             color_scheme=c("gray", "viridis"), col=NULL,
             threshold=0, swap_axes=FALSE,
             hlines=NULL, hlines_col="#B3B3B370", hlines_lwd=1, hlines_lty=1,
             vlines=NULL, vlines_col="#B3B3B370", vlines_lwd=1, vlines_lty=1,
             ...)
{
    # check inputs
    if(is.null(map)) stop("map is NULL")
    if(is.null(chr)) chr <- names(probs)[1]
    if(length(chr) > 1) {
        warning("chr should have length 1; using the first value")
        chr <- chr[1]
    }
    if(!(chr %in% names(probs))) stop("chr ", chr, " not found in probs")
    if(!(chr %in% names(map))) stop("chr ", chr, " not found in map")
    if(length(ind) > 1) {
        warning("ind should have length 1; using the first value")
        ind <- ind[1]
    }

    # pull out selected chromosomes
    probs <- probs[[chr]]
    map <- map[[chr]]
    if(dim(probs)[[3]] != length(map)) stop("Different numbers of positions in probs and map")

    # pull out individual's probs; make it positions x probs
    if(is.character(ind) && !(ind %in% rownames(probs)))
        stop("ind ", ind, " not found in probs")
    if(is.numeric(ind) && (ind < 1 || ind > nrow(probs)))
        stop("ind ", ind, " should be in the range [1, ", nrow(probs), "]")
    probs <- t(probs[ind,,])

    # pull out selected genotypes
    if(!is.null(geno)) {
        if(is.numeric(geno)) {
            if(all(geno < 0)) {
                if(any(geno > -1 | geno < -ncol(probs)))
                    stop("negative geno should be in the range [", -ncol(probs), " , -1]")
                geno <- (1:ncol(probs))[geno]
            }
            if(any(geno < 1 | geno > ncol(probs)))
                stop("numeric geno should be in the range [1, ", ncol(probs))
            geno <- colnames(probs)[geno]
        }
        if(!all(geno %in% colnames(probs)))
            stop("Not all geno in probs")
        probs <- probs[,geno,drop=FALSE]
    }

    # drop genotypes that do exceed threshold
    if(threshold > 0) {
        geno_keep <- colSums(probs >= threshold) > 0
        if(sum(geno_keep) == 0) stop("No genotype probabilities exceed the threshold")
        probs <- probs[,geno_keep,drop=FALSE]
    }

    # set up colors
    if(is.null(col)) {
        color_scheme <- match.arg(color_scheme)
        if(color_scheme=="gray") col <- gray((256:0)/256)
        else if(color_scheme=="viridis") col <- viridis_qtl2(256)
    }

    # separate positions if necessary
    tol <- 1e-6
    if(any(diff(map) < tol))
        map <- map + seq(0, tol, length.out=length(map))

    plot_genoprob_internal(probs, map, col=col, swap_axes=swap_axes,
                           hlines=hlines, hlines_col=hlines_col, hlines_lty=hlines_lty, hlines_lwd=hlines_lwd,
                           vlines=vlines, vlines_col=vlines_col, vlines_lty=vlines_lty, vlines_lwd=vlines_lwd,
                           ...)

}



# the function that does the work
plot_genoprob_internal <-
    function(probs, map, col=NULL, swap_axes=FALSE,
             zlim=c(0,1), xlab=NULL, ylab=NULL, las=NULL,
             hlines=NULL, hlines_col="gray70", hlines_lwd=1, hlines_lty=1,
             vlines=NULL, vlines_col="gray70", vlines_lwd=1, vlines_lty=1,
             mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
             ...)
{
    dots <- list(...)
    if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
    if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
    if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

    if(swap_axes) {
        probs <- t(probs)

        x <- 1:nrow(probs)
        y <- map

        if(is.null(xlab)) xlab <- ""
        if(is.null(ylab)) ylab <- "Position"

        ytick <- pretty(map)
        yticklab <- NULL

        xtick <- x
        xticklab <- rownames(probs)

        if(is.null(hlines)) hlines <- pretty(map, n=10)
        if(is.null(vlines)) vlines <- 1:nrow(probs)

        if(is.null(las)) las <- 2
    } else {
        # reverse order of genotypes, so that the first appears at the top
        probs <- probs[,ncol(probs):1,drop=FALSE]

        x <- map
        y <- 1:ncol(probs)

        if(is.null(xlab)) xlab <- "Position"
        if(is.null(ylab)) ylab <- ""

        xtick <- pretty(map)
        xticklab <- NULL

        ytick <- y
        yticklab <- colnames(probs)

        if(is.null(hlines)) hlines <- 1:nrow(probs)
        if(is.null(vlines)) vlines <- pretty(map, n=10)

        if(is.null(las)) las <- 1
    }

    image(x, y, probs,
          ylab="", yaxt="n", xlab="", xaxt="n",
          las=1, zlim=zlim, col=col, ...)

    if(dots$xaxt != "n")
        axis(side=1, at=xtick, labels=xticklab, las=las, mgp=mgp.x, tick=FALSE)
    if(dots$yaxt != "n")
        axis(side=2, at=ytick, labels=yticklab, las=las, mgp=mgp.y, tick=FALSE)


    # add grid lines
    if(!(length(hlines)==1 && is.na(hlines)))
        abline(h=hlines, lty=hlines_lty, lwd=hlines_lwd, col=hlines_col)
    if(!(length(vlines)==1 && is.na(vlines)))
        abline(v=vlines, lty=vlines_lty, lwd=vlines_lwd, col=vlines_col)

    title(xlab=xlab, mgp=mgp.x)
    title(ylab=ylab, mgp=mgp.y)
}
