#' Plot gene locations for a genomic interval
#'
#' Plot gene locations for a genomic interval, as rectangles with gene
#' symbol (and arrow indicating strand/direction) below.
#'
#' @param genes Data frame containing `start` and `stop` in
#' Mbp, `strand` (as `"-"`, `"+"`, or `NA`), and
#' `Name`.
#' @param minrow Minimum number of rows of genes in the plot
#' @param padding Proportion to pad with white space around the genes
#' @param colors Vectors of colors, used sequentially and then re-used.
#' @param scale_pos Factor by which to scale position (for example, to convert basepairs to Mbp)
#' @param start_field Character string with name of column containing the genes' start positions.
#' @param stop_field Character string with name of column containing the genes' stop positions.
#' @param strand_field Character string with name of column containing the genes' strands.
#'   (The values of the corresponding field can be character strings `"+"` or `"-"`, or numeric +1 or -1.)
#' @param name_field Character string with name of column containing the genes' names.
#' @param ... Optional arguments passed to `plot()`.
#'
#' @return None.
#'
#' @keywords hgraphics
#' @export
#' @importFrom graphics strheight strwidth text par rect abline box
#'
#' @section Hidden graphics parameters:
#' Graphics parameters can be passed via `...`. For
#' example, `xlim` to control the x-axis limits.
#' These are not included as formal
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
#'
#' # use scale_pos=1e-6 because data in bp but we want the plot in Mbp
#' plot_genes(genes, xlim=c(140, 146), scale_pos=1e-6)

# create an empty plot with test x- and y-axis limits
plot_genes <-
    function(genes, minrow=4, padding=0.2,
             colors=c("black", "red3", "green4", "blue3", "orange"),
             scale_pos=1, start_field="start", stop_field="stop",
             strand_field="strand", name_field="Name", ...)
{
    if(is.null(genes)) stop("genes is NULL")

    # make sure the columns are there
    fields <- c(start_field, stop_field, strand_field, name_field)
    fields_found <- fields %in% colnames(genes)
    if(!all(fields_found)) {
        stop("Columns not found: ", paste(fields[!fields_found], collapse=", "))
    }

    if(!is_pos_number(minrow)) stop("minrow should be a positive integer")
    if(!is_nonneg_number(padding)) stop("padding should be a non-negative number")
    if(!is_pos_number(scale_pos)) stop("scale_pos should be a positive number")

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
        genes <- genes[order(start, end),]

    # grab data
    start <- genes[,start_field]*scale_pos # convert to Mbp
    end <- genes[,stop_field]*scale_pos    # convert to Mbp
    strand <- as.character(genes[,strand_field])
    name <- as.character(genes[,name_field])

    # deal with +/- 1 for strand
    strand[!is.na(strand) & (strand=="1" | strand=="+1")] <- "+"
    strand[!is.na(strand) & strand=="-1"] <- "-"

    plot_genes_internal <-
        function(xlab=NULL, xaxs="i",
                 bgcolor="gray92", xat=NULL,
                 mgp=c(1.6,0.2,0), xlim=NULL,
                 vlines=NULL, vlines_col="white",
                 vlines_lwd=1, vlines_lty=1,
                 xaxt="s", ylab="", ...)
        {
            if(is.null(xlab)) {
                if(length(unique(genes$chr)) == 1)
                    xlab <- paste("Chr", genes$chr[1], "position (Mbp)")
                else xlab <- "Position (Mbp)"
            }

            if(is.null(xlim)) {
                xlim <- range(c(start, end), na.rm=TRUE)
            }

            plot(0, 0, type="n",
                 xlim=xlim,   xlab=xlab, xaxs=xaxs, xaxt="n",
                 ylim=c(1,0), ylab=ylab, yaxs="i", yaxt="n", mgp=mgp, ...)

            # gray background
            u <- par("usr")
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border="black")

            # axis
            if(is.null(xat)) xat <- pretty(xlim)
            if((length(xat) > 1 || !is.na(xat)) && xaxt != "n")
                axis(side=1, at=xat, mgp=mgp, tick=FALSE)

            # vertical lines
            if(is.null(vlines)) vlines <- xat
            if(length(vlines) > 1 || !is.na(vlines))
                abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)

            box()
        }

    plot_genes_internal(...)

    # drop genes that are not in plotting region
    u <- par("usr")
    omit <- (end < u[1] | start > u[2])
    if(any(omit)) {
        keep <- !omit
        start <- start[keep]
        end <- end[keep]
        strand <- strand[keep]
        name <- name[keep]
    }
    if(length(start) == 0) # no genes to plot
        return(invisible(NULL))

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
