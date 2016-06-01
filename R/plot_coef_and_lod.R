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
    # also, match markers and use map in coefficients object
    mar_in_coef <- rownames(x$coef)
    mar_in_scan1 <- rownames(scan1_output$lod)
    scan1_output$lod <- scan1_output$lod[mar_in_scan1 %in% mar_in_coef, , drop=FALSE]
    scan1_output$map <- list("1"=x$map)
    mis_mar <- !(mar_in_coef %in% mar_in_scan1)
    if(any(mis_mar)) {
        n_new <- sum(mis_mar)
        new_lod <- matrix(NA, nrow=sum(mis_mar), ncol=scan1_output$lod)
        rownames(new_lod) <- mar_in_coef[mis_mar]
        scan1_output$lod <- rbind(scan1_output$lod, new_lod)[mar_in_coef,]
    }


    # 2 x 1 panels; adjust margins
    old_mfrow <- par("mfrow")
    old_mar <- par("mar")
    on.exit({ cat("exiting\n"); par(mfrow=old_mfrow, mar=old_mar)})
    par(mfrow=c(2,1))
    top_mar <- bottom_mar <- old_mar
    top_mar[1] <- 0.1
    bottom_mar[3] <- 0.1

    par(mar=top_mar)
    plot_coef(x, columns, col, scan1_output=NULL,
              add=FALSE, gap=gap, ylim=ylim, bgcolor=bgcolor,
              altbgcolor=altbgcolor, ylab=ylab,
              xaxt="n", vlines=vlines, vlines.col=vlines.col,
              vlines.lwd=vlines.lwd, vlines.lty=vlines.lty, ...)

    par(mar=bottom_mar)
    plot_scan1(scan1_output, lodcolumn=1,
               add=FALSE, gap=gap, vlines=vlines, vlines.col=vlines.col,
               vlines.lwd=vlines.lwd, vlines.lty=vlines.lty, ...)
}
