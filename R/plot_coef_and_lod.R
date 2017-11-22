# two-panel plot with both coefficients and LOD scores
# (for a single chromosome)
#
# calls plot_coef and plot_scan1
# internal function that is called by plot_coef
plot_coef_and_lod <-
    function(x, map, columns=NULL, col=NULL, scan1_output,
             gap=25, ylim=NULL, bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects",
             ylab_lod="LOD score", ylim_lod=NULL, col_lod="slateblue",
             xaxt=NULL,
             vlines=NULL, vlines_col="white", vlines_lwd=1, vlines_lty=1,
             top_panel_prop=0.65, ...)
{
    if(is.null(map)) stop("map is NULL")

    # align scan1 output and map
    tmp <- align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map

    if(nrow(x) != length(unlist(map)))
        stop("nrow(x) [", nrow(x), "] != number of positions in map [",
             length(unlist(map)), "]")

    # also, match markers and use map in coefficients object
    mar_in_coef <- rownames(x)
    mar_in_scan1 <- rownames(scan1_output)
    scan1_output <- unclass(scan1_output)[mar_in_scan1 %in% mar_in_coef, , drop=FALSE]
    mis_mar <- !(mar_in_coef %in% mar_in_scan1)
    if(any(mis_mar)) {
        n_new <- sum(mis_mar)
        new_lod <- matrix(NA, nrow=sum(mis_mar), ncol=scan1_output$lod)
        rownames(new_lod) <- mar_in_coef[mis_mar]
        scan1_output <- rbind(scan1_output, new_lod)[mar_in_coef,]
    }

    # 2 x 1 panels; adjust margins
    old_mfrow <- par("mfrow")
    old_mar <- par("mar")
    on.exit(par(mfrow=old_mfrow, mar=old_mar))
    layout(rbind(1,2), heights=c(top_panel_prop, 1-top_panel_prop))
    top_mar <- bottom_mar <- old_mar
    top_mar[1] <- 0.1
    bottom_mar[3] <- 0.1

    par(mar=top_mar)
    plot_coef(x, map, columns=columns, col=col, scan1_output=NULL,
              add=FALSE, gap=gap, ylim=ylim, bgcolor=bgcolor,
              altbgcolor=altbgcolor, ylab=ylab,
              xaxt="n", vlines=vlines, vlines_col=vlines_col,
              vlines_lwd=vlines_lwd, vlines_lty=vlines_lty, ...)

    par(mar=bottom_mar)
    plot_scan1(scan1_output, map, lodcolumn=1, col=col_lod, ylab=ylab_lod,
               add=FALSE, gap=gap, bgcolor=bgcolor, altbgcolor=altbgcolor,
               vlines=vlines, vlines_col=vlines_col,
               vlines_lwd=vlines_lwd, vlines_lty=vlines_lty, ...)
}
