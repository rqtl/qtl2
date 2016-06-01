# two-panel plot with both coefficients and LOD scores
# (for a single chromosome)
#
# calls plot_coef and plot_scan1
# internal function that is called by plot_coef
plot_coef_and_lod <-
    function(x, columns=NULL, col=NULL, scan1_output,
             gap=25, ylim=NULL, bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects", ylab_lod="LOD score", ylim_lod=NULL,
             xlab="Chromosome", xaxt=NULL,
             vlines=NULL, vlines.col="white", vlines.lwd=1, vlines.lty=1,
             ...)
{
    # assume one lod score
    # if multiple, take first

    # use same chromosomes as in coefficients
    # also, match markers and use map in coefficients object
    # (NAs for missing markers?)

    # 2 x 1 panels
    old_mfrow <- par("mfrow")
    on.exit(par(mfrow=old_mfrow))
    par(mfrow=c(2,1))

    # adjust margins
    old_mar <- par("mar")
    on.exit(par(mar=old_mar))
    top_mar <- bottom_mar <- old_mar
    top_mar[1] <- 0
    bottom_mar[3] <- 0

    par(mar=top_mar)
    plot_coef(x, columns=columns, col=col, scan1_output=NULL,
              add=FALSE, gap=gap, ylim=ylim, bgcolor=bgcolor,
              altbgcolor=altbgcolor, ylab=ylab,
              xaxt="n", vlines=vlines, col=vlines.col, lwd=vlines.lwd,
              lty=vlines.lty, ...)

    par(mar=bottom_mar)
    plot_scan1(scan1_output, lodcolumn=1, chr=NULL, add=FALSE,
               gap=gap, vlines=vlines, vlines.col=vlines.col,
               vlines.lwd=vlines.lwd, vlines.lty=vlines.lty, ...)

}
