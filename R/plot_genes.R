#' Plot gene locations for a genomic interval
#'
#' Plot gene locations for a genomic interval, as rectangles with gene
#' symbol (and arrow indicating strand/direction) below.
#'
#' @param genes Data frame containing \code{start} and \code{stop} in
#' bp, \code{strand} (as \code{"-"}, \code{"+"}, or \code{NA}), and
#' \code{Name}.
#' @param xlim x-axis limits (in Mbp)
#' @param minrow Minimum number of rows of genes in the plot
#' @param padding Proportion to pad with white space around the genes
#' @param colors Vectors of colors, used sequentially and then re-used.
#' @param scale_pos Factor by which to scale position (default converts bp to Mbp)
#' @param start_field Character string with name of column containing the genes' start positions.
#' @param stop_field Character string with name of column containing the genes' stop positions.
#' @param strand_field Character string with name of column containing the genes' strands.
#' @param name_field Character string with name of column containing the genes' names.
#' @param ... Optional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return None.
#'
#' @keywords hgraphics
#' @export
#' @importFrom graphics strheight strwidth text plot par rect abline box
#'
#' @examples
#' genes <- data.frame(chr = c("6", "6", "6", "6", "6", "6", "6", "6"),
#'                     start = c(139988753, 140680185, 141708118, 142234227, 142587862,
#'                               143232344, 144398099, 144993835),
#'                     stop  = c(140041457, 140826797, 141773810, 142322981, 142702315,
#'                               143260627, 144399821, 145076184),
#'                     strand = c("-", "+", "-", "-", "-", NA, "+", "-"),
#'                     Name = c("Plcz1", "Gm30215", "Gm5724", "Slco1a5", "Abcc9",
#'                              "4930407I02Rik", "Gm31777", "Bcat1"),
#'                     stringsAsFactors=FALSE)
#' plot_genes(genes, xlim=c(140, 146))

# create an empty plot with test x- and y-axis limits
plot_genes <-
    function(genes, xlim=NULL, minrow=4, padding=0.2,
             colors=c("black", "red3", "green4", "blue3", "orange"),
             scale_pos=1e-6, start_field="start", stop_field="stop",
             strand_field="strand", name_field="Name", ...)
{
    # make sure the columns are there
    fields <- c(start_field, stop_field, strand_field, name_field)
    fields_found <- fields %in% colnames(genes)
    if(!all(fields_found)) {
        stop("Columns not found: ", paste(fields[!fields_found], collapse=", "))
    }

    # grab just the start and stop
    start <- genes[,start_field]
    end <- genes[,stop_field]

    # drop genes with missing start or stop
    missing_pos <- is.na(start) | is.na(end)
    if(any(missing_pos)) {
        warning("Dropping ", sum(missing_pos), " rows with missing positions")
        genes <- genes[!missing_pos, , drop=FALSE]
    }

    # make sure genes are ordered by their start values
    if(any(diff(start) < 0))
        genes <- genes[order(start, stop),]

    # grab data
    start <- genes[,start_field]*scale_pos # convert to Mbp
    end <- genes[,stop_field]*scale_pos    # convert to Mbp
    strand <- as.character(genes[,strand_field])
    name <- as.character(genes[,name_field])

    if(is.null(xlim)) {
        xlim <- range(c(start, end), na.rm=TRUE)
    }

    internal_plot_genes <-
        function(xlab="Position (Mbp)", xaxs="i",
                 bgcolor="gray92", xat=NULL,
                 mgp=c(1.6,0.2,0),
                 vlines=NULL, vlines_col="white",
                 vlines_lwd=1, vlines_lty=1)
        {
            plot(0, 0, type="n",
                 xlim=xlim,   xlab=xlab, xaxs=xaxs, xaxt="n",
                 ylim=c(1,0), ylab="", yaxs="i", yaxt="n", mgp=mgp)

            # gray background
            u <- par("usr")
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border="black")

            # axis
            if(is.null(xat)) xat <- pretty(xlim)
            if(length(xat) > 1 || !is.na(xat))
                axis(side=1, at=xat, mgp=mgp, tick=FALSE)

            # vertical lines
            if(is.null(vlines)) vlines <- xat
            if(length(vlines) > 1 || !is.na(vlines))
                abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)

            box()
        }

    internal_plot_genes(...)

    # missing names: use ?
    name[is.na(name)] <- "?"

    # arrow annotation re direction, to place after gene name
    dir_symbol <- rep(' ', length(name))
    right <- !is.na(strand) & strand == "+"
    if(any(right))
        dir_symbol[right] <- expression(phantom('') %->% phantom(''))
    left <- !is.na(strand) & strand == "-"
    if(any(left))
        dir_symbol[left] <- expression(phantom('') %<-% phantom(''))

    # initial determination of text size
    maxy <- minrow
    height <- 1/maxy
    text_cex <- 1

    # adjust text size and determine vertical location of genes
    for(it in 1:2) { # go through all of this twice

        while(max(abs(strheight(name, cex=text_cex))) > height*(1-padding)) {
            text_cex <- text_cex * 0.99
        }

        # horizontal padding
        space <- strwidth(' ', cex=text_cex)

        # figure out how to arrange genes vertically
        #     + number of rows of genes
        # (function defined in src/arrange_genes.cpp)
        y <- arrange_genes(start, end + space + strwidth(name, cex=text_cex) + strwidth(dir_symbol, cex=text_cex))

        maxy <- max(c(y, minrow))
        height <- 1/maxy
    }

    ypos <- seq(height/2, by=height, length=maxy)
    y <- ypos[y]
    rect_height <- height*(1-padding)
    rect_top <- y - rect_height/2
    rect_bottom <- y + rect_height/2

    colors <- rep(colors, length(y))
    for(i in seq(along=start)) {
        rect(start[i], rect_top[i],
             end[i],   rect_bottom[i],
             col=colors[i], border=colors[i],
             lend=1, ljoin=1)
        text(end[i] + space, y[i],
             name[i], adj=c(0, 0.5), col=colors[i],
             cex=text_cex)
        if(!is.na(strand[i]) && (strand[i] == "+" || strand[i] == '-'))
            text(end[i] + space + strwidth(name[i], cex=text_cex), y[i],
                 dir_symbol[i], adj=c(0, 0.5), col=colors[i],
                 cex=text_cex)
    }
}
