#' plot strain distribution patterns for SNPs
#'
#' plot the strain distribution patterns of SNPs
#' using tracks of tick-marks for each founder strain
#'
#' @param pos vector of SNP positions
#' @param sdp vector of strain distribution patterns (as integers)
#' @param labels names of the strains
#' @param ... additional graphic arguments
#'
#' @return None.
#'
#' @keywords hgraphics
#' @export
#' @importFrom graphics segments axis title rect box abline
#'
#' @seealso [calc_sdp()], [invert_sdp()]
#'
#' @details
#' Additional arguments, such as `xlab`, `ylab`, `xlim`, and `main`,
#' are passed via `...`; also `bgcolor` to control the color of the
#' background, and `col` and `lwd` to control the color and thickness
#' of the tick marks.
#'
#' @examples
#' n_tick <- 50
#' plot_sdp(runif(n_tick, 0, 100), sample(0:255, n_tick, replace=TRUE))

plot_sdp <-
    function(pos, sdp, labels=names(qtl2::CCcolors), ...)
{
    n_str <- length(labels)

    stopifnot(length(pos) == length(sdp))
    stopifnot(all(sdp < 2^n_str & sdp >= 0))

    alleles <- invert_sdp(sdp, n_str)

    y <- seq_len(n_str)

    plot_sdp_internal <-
        function(xlim=range(pos), ylim=c(max(y)+0.5, min(y)-0.5),
                 xlab="Position (Mbp)", ylab="",
                 xaxs="i", yaxs="i",
                 main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
                 mgp=NULL, las=1,
                 hlines=NULL, hlines_col="black", hlines_lwd=1, hlines_lty=1,
                 vlines=NULL, vlines_col="white", vlines_lwd=1, vlines_lty=1,
                 bgcolor="gray90", lwd=2, col="darkslateblue", sub="", ...)
        {
            dots <- list(...)

            # margin parameters
            if(!is.null(mgp)) mgp.x <- mgp.y <- mgp

            # make basic plot
            plot(0, 0, xlab="", ylab="", xlim=xlim, ylim=ylim,
                 xaxs=xaxs, yaxs=yaxs, xaxt="n", yaxt="n", type="n",
                 main=main, sub=sub)

            # add background rectangle
            u <- par("usr")
            if(!is.null(bgcolor))
                rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)

            # include axis labels?
            if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
            if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

            if(dots$xaxt != "n") {
                axis(side=1, at=pretty(xlim), mgp=mgp.x, las=las, tick=FALSE)
            }

            if(dots$yaxt != "n") {
                axis(side=2, at=y, labels, mgp=mgp.y, las=las, tick=FALSE)
            }

            # add axis titles
            title(xlab=xlab, mgp=mgp.x)
            title(ylab=ylab, mgp=mgp.y)

            # grid lines
            if(!(length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
                if(is.null(hlines)) hlines <- seq(0.5, n_str+0.5)
                abline(h=hlines, col=hlines_col, lwd=hlines_lwd, lty=hlines_lty)
            }
            if(!(length(vlines)==1 && is.na(vlines))) { # if vlines==NA (or mult chr), skip lines
                if(is.null(vlines)) vlines <- pretty(xlim)
                abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
            }


            for(i in seq_len(n_str)) {
                if(all(alleles[,i] == 1)) next
                wh <- which(alleles[,i]==3)
                yy <- rep(y[i], length(wh))
                segments(pos[wh], yy-0.5, pos[wh], yy+0.5, lwd=lwd, col=col)
            }

            box()
        }

    plot_sdp_internal(...)

}
