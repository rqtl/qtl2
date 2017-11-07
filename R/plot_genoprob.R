#' Plot genotype probabilities for one individual on one chromosome.
#'
#' Plot the genotype probabilities for one individual on one chromosome, as a heat map.
#'
#' @param probs Genotype probabilities (as produced by \code{\link[qtl2geno]{calc_genoprob}})
#' or allele dosages (as produced by \code{\link[qtl2geno]{genoprob_to_alleleprob}}).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosome to plot; a single character string.
#' @param geno Optional vector of genotypes or alleles to be shown.
#' @param colscheme Color scheme for the heatmap (ignored if \code{col} is provided).
#' @param col Optional vector of colors for the heatmap.
#' @param threshold Threshold for genotype probabilities; only genotypes that achieve
#' this value somewhere on the chromosome will be shown.
#' @param transpose If TRUE, swap the axes, so that the genotypes are
#' on the x-axis and the chromosome position is on the y-axis.
#' @param hlines Position of horizontal grid lines (use \code{NA} to avoid lines).
#' @param hlines.col Color of horizontal grid lines.
#' @param hlines.lwd Line width of horizontal grid lines.
#' @param hlines.lty Line type of horizontal grid lines.
#' @param vlines Position of vertical grid lines (use \code{NA} to avoid lines).
#' @param vlines.col Color of vertical grid lines.
#' @param vlines.lwd Line width of vertical grid lines.
#' @param vlines.lty Line type of vertical grid lines.
#'
#' @param ... Additional graphics parameters passed to \code{\link[graphics]{image}}.
#'
#' @export
#' @importFrom graphics image par axis title box
plot_genoprob <-
    function(probs, map, ind=1, chr=NULL, geno=NULL,
             color_scheme=c("gray", "viridis"), col=NULL,
             threshold=0, transpose=FALSE,
             hlines=NULL, hlines.col="#B3B3B370", hlines.lwd=1, hlines.lty=1,
             vlines=NULL, vlines.col="#B3B3B370", vlines.lwd=1, vlines.lty=1,
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
    # reverse order of genotypes, so that the first appears at the top
    probs <- probs[,ncol(probs):1,drop=FALSE]

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

    # the function that does the work
    plot_genoprob_internal <-
        function(probs, map, col=NULL, transpose=FALSE,
                 zlim=c(0,1), xlab=NULL, ylab=NULL, main="", las=NULL,
                 hlines=NULL, hlines.col="gray70", hlines.lwd=1, hlines.lty=1,
                 vlines=NULL, vlines.col="gray70", vlines.lwd=1, vlines.lty=1,
                 mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
                 ...)
    {
        dots <- list(...)
        if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
        if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
        if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

        if(transpose) {
            # reverse order of genotypes back so that first appears at left
            #    then transpose
            probs <- t(probs[,ncol(probs):1,drop=FALSE])

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
            abline(h=hlines, lty=hlines.lty, lwd=hlines.lwd, col=hlines.col)
        if(!(length(vlines)==1 && is.na(vlines)))
            abline(v=vlines, lty=vlines.lty, lwd=vlines.lwd, col=vlines.col)

        title(xlab=xlab, mgp=mgp.x)
        title(ylab=ylab, mgp=mgp.y)
    }

    plot_genoprob_internal(probs, map, col=col, transpose=transpose,
                           hlines=hlines, hlines.col=hlines.col, hlines.lty=hlines.lty, hlines.lwd=hlines.lwd,
                           vlines=vlines, vlines.col=vlines.col, vlines.lty=vlines.lty, vlines.lwd=vlines.lwd,
                           ...)

}
