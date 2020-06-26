#' Plot LOD scores vs QTL peak locations
#'
#' Create a scatterplot of LOD scores vs QTL peak locations (possibly with intervals) for multiple traits.
#'
#' @param peaks Data frame such as that produced by
#'     [find_peaks()]) containing columns
#'     `chr`, `pos`, `lodindex`, and `lodcolumn`.
#'     May also contain columns `ci_lo` and `ci_hi`, in
#'     which case intervals will be plotted.
#' @param map Marker map, used to get chromosome lengths (and start
#'     and end positions).
#' @param chr Selected chromosomes to plot; a vector of character
#'     strings.
#' @param gap Gap between chromosomes. The default is 1% of the total genome length.
#' @param intervals If TRUE and `peaks` contains QTL intervals, plot the intervals.
#' @param ... Additional graphics parameters
#'
#' @seealso [find_peaks()], [plot_peaks()]
#'
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color and
#' `altbgcolor` to control the background color on alternate chromosomes.
#' These are not included as formal parameters in order to avoid
#' cluttering the function definition.
#'
#' @export
#' @importFrom graphics segments abline par axis title box rect
#'
#' @return None.
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c(1,2,7,8,9,13,15,16,19)]}
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
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # find peaks above lod=3.5 (and calculate 1.5-LOD support intervals)
#' peaks <- find_peaks(out, map, threshold=3.5, drop=1.5)
#'
#' plot_lodpeaks(peaks, map)

plot_lodpeaks <-
    function(peaks, map, chr=NULL, gap=NULL, intervals=FALSE, ...)
{
    if(is.null(peaks)) stop("peaks is NULL")
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

    plot_lodpeaks_internal <-
        function(peaks, map, gap, intervals=FALSE,
                 bgcolor="gray90", altbgcolor="gray85",
                 lwd=2, pch=21, col="black", bg="slateblue",
                 xlab=NULL, ylab="LOD score",
                 xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
                 main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
                 mgp=NULL, las=1,
                 hlines=NULL, hlines_col="white", hlines_lwd=1, hlines_lty=1,
                 vlines=NULL, vlines_col="white", vlines_lwd=1, vlines_lty=1,
                 lend=1, ljoin=1, sub="", ...)
        {
            dots <- list(...)
            onechr <- (length(map)==1) # single chromosome

            if(is.null(ylim)) {
                ylim <- c(0, max(peaks$lod, na.rm=TRUE)*1.05)
            }

            # plot empty frame
            emptylod <- as.matrix(rep(NA, sum(sapply(map, length))))
            rownames(emptylod) <- unlist(lapply(map, names))
            plot_scan1(emptylod, map, gap=gap, type="n",
                       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs,
                       bgcolor=bgcolor, altbgcolor=altbgcolor,
                       main=main, sub=sub, las=las,
                       hlines=hlines, hlines_col=hlines_col, hlines_lwd=hlines_lwd, hlines_lty=hlines_lty,
                       vlines=vlines, vlines_col=vlines_col, vlines_lwd=vlines_lwd, vlines_lty=vlines_lty)

            # add QTL intervals
            if(intervals && all(c("ci_lo", "ci_hi") %in% names(peaks))) {
                xpos_lo <- xpos_scan1(map, gap=gap, thechr=peaks$chr, thepos=peaks$ci_lo)
                xpos_hi <- xpos_scan1(map, gap=gap, thechr=peaks$chr, thepos=peaks$ci_hi)

                segments(xpos_lo, peaks$lod, xpos_hi, peaks$lod, lend=lend, ljoin=ljoin)
            }

            # add points at LOD scores
            xpos <- xpos_scan1(map, gap=gap, thechr=peaks$chr, thepos=peaks$pos)
            points(xpos, peaks$lod, pch=pch, col=col, bg=bg, ...)
        }

    plot_lodpeaks_internal(peaks=peaks, map=map, gap=gap, intervals=intervals, ...)

}
