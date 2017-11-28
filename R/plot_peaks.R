#' Plot QTL peak locations
#'
#' Plot QTL peak locations (possibly with intervals) for multiple traits.
#'
#' @md
#'
#' @param peaks Data frame such as that produced by
#'     [qtl2scan::find_peaks()]) containing columns
#'     `chr`, `pos`, `lodindex`, and `lodcolumn`.
#'     May also contain columns `ci_lo` and `ci_hi`, in
#'     which case intervals will be plotted.
#' @param map Marker map, used to get chromosome lengths (and start
#'     and end positions).
#' @param chr Selected chromosomes to plot; a vector of character
#'     strings.
#' @param tick_height Height of tick marks at the peaks (a number between 0 and 1).
#' @param gap Gap between chromosomes.
#' @param ... Additional graphics parameters
#'
#' @seealso [qtl2scan::find_peaks()]
#'
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color and
#' `altbgcolor` to control the background color on alternate chromosomes.
#' These are not included as formal parameters in order to avoid
#' cluttering the function definition.
#'
#' @export
#' @importFrom graphics plot segments abline par axis title box rect
#'
#' @return None.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
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
#' library(qtl2scan)
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # find peaks above lod=3.5 (and calculate 1.5-LOD support intervals)
#' peaks <- find_peaks(out, map, threshold=3.5, drop=1.5)
#'
#' plot_peaks(peaks, map)

plot_peaks <-
    function(peaks, map, chr=NULL, tick_height = 0.3,
             gap=25, ...)
{
    if(is.null(map)) stop("map is NULL")

    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        map <- map[chri]
        peaks <- peaks[peaks$chr %in% names(map),,drop=FALSE]
    }

    if(nrow(peaks) == 0)
        stop("There are no peaks on the selected chromosomes")

    plot_peaks_internal <-
        function(peaks, map, tick_height,
                 gap, bgcolor="gray90", altbgcolor="gray85",
                 lwd=2, col="slateblue", xlab=NULL, ylab="",
                 xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
                 main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
                 mgp=NULL, las=1, lend=1, ljoin=1,
                 hlines=NULL, hlines_col="white", hlines_lwd=1, hlines_lty=1,
                 vlines=NULL, vlines_col="white", vlines_lwd=1, vlines_lty=1,
                 ...)
        {
            dots <- list(...)
            onechr <- (length(map)==1) # single chromosome

            xpos <- map_to_xpos(map, gap)
            chrbound <- map_to_boundaries(map, gap)

            unique_lodindex <- sort(unique(peaks$lodindex))

            if(is.null(ylim))
                ylim <- c(length(unique_lodindex)+0.5, 0.5)

            if(is.null(xlim)) {
                xlim <- range(xpos, na.rm=TRUE)
                if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
            }

            if(is.null(xlab)) {
                if(onechr) {
                    if(names(map) == " ") xlab <- "Position"
                    else xlab <- paste("Chr", names(map), "position")
                }
                else xlab <- "Chromosome"
            }

            # margin parameters
            if(!is.null(mgp)) mgp.x <- mgp.y <- mgp

            # make basic plot
            plot(0,0, xlab="", ylab="", xlim=xlim, ylim=ylim,
                 xaxs=xaxs, yaxs=yaxs, xaxt="n", yaxt="n", type="n",
                 main=main)

            # add background rectangles
            u <- par("usr")
            if(!is.null(bgcolor))
                rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
            if(!is.null(altbgcolor) && !onechr) {
                for(i in seq(2, ncol(chrbound), by=2))
                    rect(chrbound[1,i], u[3], chrbound[2,i], u[4], col=altbgcolor, border=NA)
            }

            # include axis labels?
            if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
            if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

            # add x axis unless par(xaxt="n")
            if(dots$xaxt != "n") {
                if(onechr) {
                    axis(side=1, at=pretty(xlim), mgp=mgp.x, las=las, tick=FALSE)
                }
                else {
                    loc <- colMeans(chrbound)
                    odd <- seq(1, length(map), by=2)
                    even <- seq(2, length(map), by=2)
                    axis(side=1, at=loc[odd], names(map)[odd],
                         mgp=mgp.x, las=las, tick=FALSE)
                    axis(side=1, at=loc[even], names(map)[even],
                         mgp=mgp.x, las=las, tick=FALSE)
                }
            }

            # add y axis unless par(yaxt="n")
            if(dots$yaxt != "n") {
                axis(side=2, at=seq(along=unique_lodindex),
                     peaks$lodcolumn[match(unique_lodindex, peaks$lodindex)], # <- the y-axis labels
                     mgp=mgp.y, las=las, tick=FALSE)
            }
            dots$xaxt <- dots$yaxt <- NULL # delete those

            # x and y axis labels
            title(xlab=xlab, mgp=mgp.x)
            title(ylab=ylab, mgp=mgp.y)

            # grid lines
            if(onechr && !(length(vlines)==1 && is.na(vlines))) { # if vlines==NA (or mult chr), skip lines
                if(is.null(vlines)) vlines <- pretty(xlim)
                abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
            }
            if(!(length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
                if(is.null(hlines)) hlines <- seq_along(unique_lodindex)
                abline(h=hlines, col=hlines_col, lwd=hlines_lwd, lty=hlines_lty)
            }

            y <- match(peaks$lodindex, unique_lodindex)
            x <- xpos_scan1(map, names(map), gap, peaks$chr, peaks$pos)

            segments(x, y-tick_height/2, x, y+tick_height/2,
                     col=col, lwd=lwd, lend=lend, ljoin=ljoin)

            if("ci_lo" %in% names(peaks) && "ci_hi" %in% names(peaks)) {
                x_lo <- xpos_scan1(map, names(map), gap, peaks$chr, peaks$ci_lo)
                x_hi <- xpos_scan1(map, names(map), gap, peaks$chr, peaks$ci_hi)

                segments(x_lo, y, x_hi, y, col=col, lwd=lwd,
                         lend=lend, ljoin=ljoin)
            }

        }

    plot_peaks_internal(peaks=peaks, map=map, tick_height=tick_height,
                        gap=gap, ...)

}
