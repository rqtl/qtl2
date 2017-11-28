# two-panel plot with both snp asso and genes
# (for a single chromosome)
#
# calls plot_snpasso and plot_genes
# internal function that is called by plot_snpasso
plot_snpasso_and_genes <-
    function(scan1output, snpinfo, show_all_snps=TRUE,
             drop_hilit=NA, col_hilit="violetred", col="darkslateblue",
             gap=25, minlod=0,
             genes, minrow=4, padding=0.2,
             colors=c("black", "red3", "green4", "blue3", "orange"),
             scale_pos=1, start_field="start", stop_field="stop",
             strand_field="strand", name_field="Name",
             top_panel_prop=0.65, xlim=NULL, xaxt="s",
             xlab=NULL, ...)
{
    # 2 x 1 panels; adjust margins
    old_mfrow <- par("mfrow")
    old_mar <- par("mar")
    on.exit(par(mfrow=old_mfrow, mar=old_mar))
    layout(rbind(1,2), heights=c(top_panel_prop, 1-top_panel_prop))
    top_mar <- bottom_mar <- old_mar
    top_mar[1] <- 0.1
    bottom_mar[3] <- 0.1

    if(is.null(xlim)) xlim <- range(snpinfo$pos)

    if(is.null(xlab)) {
        if(length(unique(snpinfo$chr))==1)
            xlab <- paste("Chr", snpinfo$chr[1], "position (Mbp)")
        else
            xlab <- "Position (Mbp)"
    }

    par(mar=top_mar)
    plot_snpasso(scan1output, snpinfo, show_all_snps=show_all_snps,
                 drop_hilit=drop_hilit, col_hilit=col_hilit, col=col,
                 gap=gap, minlod=minlod, xlim=xlim, xaxt="n", xlab="",
                 ...)

    par(mar=bottom_mar)
    plot_genes(genes, minrow=minrow, padding=padding, colors=colors,
               scale_pos=scale_pos, start_field=start_field, stop_field=stop_field,
               strand_field=strand_field, name_field=name_field, xlim=xlim,
               xaxt=xaxt, xlab=xlab, ...)

}
