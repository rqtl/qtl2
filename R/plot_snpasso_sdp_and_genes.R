# three-panel plot with both snp asso, SDP, and genes
# (for a single chromosome)
#
# calls plot_snpasso, plot_sdp, and plot_genes
# internal function that is called by plot_snpasso
plot_snpasso_sdp_genes <-
    function(scan1output, snpinfo, show_all_snps=TRUE,
             drop_hilit=NA, col_hilit="violetred", col="darkslateblue",
             gap=NULL, minlod=0,
             genes, minrow=4, padding=0.2,
             colors=c("black", "red3", "green4", "blue3", "orange"),
             scale_pos=1, start_field="start", stop_field="stop",
             strand_field="strand", name_field="Name",
             panel_prop=c(0.2,0.45,0.35), xlim=NULL, xaxt="s",
             xlab=NULL, main="", sub="",
             strain_labels=names(qtl2::CCcolors), ...)
{
    # 3 x 1 panels; adjust margins
    old_mfrow <- par("mfrow")
    old_mar <- par("mar")
    on.exit(par(mfrow=old_mfrow, mar=old_mar))
    layout(rbind(1,2,3), heights=panel_prop)
    top_mar <- middle_mar <- bottom_mar <- old_mar
    top_mar[1] <- middle_mar[1] <- middle_mar[3] <- bottom_mar[3] <- 0.1

    if(is.null(xlim)) xlim <- range(snpinfo$pos)

    if(is.null(xlab)) {
        if(length(unique(snpinfo$chr))==1)
            xlab <- paste("Chr", snpinfo$chr[1], "position (Mbp)")
        else
            xlab <- "Position (Mbp)"
    }


    # determine snps to show in SDP plot
    # maybe expand snp info
    map <- snpinfo_to_map(snpinfo)
    if(show_all_snps) {
        tmp <- expand_snp_results(scan1output, map, snpinfo)
        scan1output <- tmp$lod
        map <- tmp$map
    }
    if(is.na(drop_hilit)) drop_hilit <- Inf
    snps2show <- rownames(scan1output)[max(scan1output[,1]) - scan1output[,1] <= drop_hilit]
    snpinfo_sub <- snpinfo[snpinfo$snp %in% snps2show,,drop=FALSE]

    par(mar=top_mar)
    plot_sdp(snpinfo_sub$pos, snpinfo_sub$sdp, strain_labels=strain_labels,
             xlim=xlim, xaxt="n", xlab="", main=main, ...)

    par(mar=middle_mar)
    plot_snpasso(scan1output, snpinfo, show_all_snps=show_all_snps,
                 drop_hilit=drop_hilit, col_hilit=col_hilit, col=col,
                 gap=gap, minlod=minlod, xlim=xlim, xaxt="n", xlab="",
                 ...)

    par(mar=bottom_mar)
    plot_genes(genes, minrow=minrow, padding=padding, colors=colors,
               scale_pos=scale_pos, start_field=start_field, stop_field=stop_field,
               strand_field=strand_field, name_field=name_field, xlim=xlim,
               xaxt=xaxt, xlab=xlab, sub=sub, ...)

}
