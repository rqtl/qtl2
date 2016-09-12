#' Plot gene locations for a genomic interval
#'
#' Plot gene locations for a genomic interval, as rectangles with gene
#' symbol (and arrow indicating strand/direction) below.
#'
#' @param genes Data frame containing \code{start} and \code{stop} in
#' bp, \code{strand} (as \code{"-"}, \code{"+"}, or \code{NA}), and
#' \code{Name}.
#' @param xlim x-axis limits (in Mbp)
#' @param minrow Minimum number of rows of genes
#' @param padding Proportion to pad with white space around the genes
#' @param colors Vectors of colors, used sequentially and then re-used.
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
             ...)
{
    # grab data
    start <- genes$start/10^6 # convert to Mbp
    end <- genes$stop/10^6   # convert to Mbp
    strand <- as.character(genes$strand)
    name <- as.character(genes$Name)

    if(is.null(xlim)) {
        xlim <- range(c(start, end), na.rm=TRUE)
    }

    internal_plot_genes <-
        function(xlab="Position (Mbp)", xaxs="i",
                 bgcolor="gray92", xat=NULL,
                 mgp=c(0,0.2,0),
                 vlines=NULL, vlines.col="white",
                 vlines.lwd=1, vlines.lty=1)
        {
            plot(0, 0, type="n",
                 xlim=xlim,   xlab=xlab, xaxs=xaxs, xaxt="n",
                 ylim=c(1,0), ylab="", yaxs="i", yaxt="n")

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
                abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)

            box()
        }

    internal_plot_genes(...)

    # missing names: use ?
    name[is.na(name)] <- "?"

    # right and left arrows (unicode symbol codes)
    right_arrow <- "\u2192"
    left_arrow <- "\u2190"

    # add direction to gene name
    right <- !is.na(strand) & strand == "+"
    left <- !is.na(strand) & strand == "-"
    if(any(right))
        name[right] <- paste(name[right], right_arrow)
    if(any(left))
        name[left] <- paste(name[left], left_arrow)

    # initial determination of text size
    maxy <- minrow
    height <- 1/maxy
    text_cex <- 1

    # adjust text size and determine vertical location of genes
    for(it in 1:2) { # go through all of this twice

        while(max(abs(strheight(name, cex=text_cex))) > height/2*(1-padding)) {
            text_cex <- text_cex * 0.99
        }

        # horizontal padding
        xpad <- strwidth("m", cex=text_cex)

        # figure out how to arrange genes vertically
        #   + number of rows of genes
        n <- nrow(genes)
        y <- rep(NA, n)
        maxy <- y[1] <- 1
        maxx <- max(c(end[1], start[1] + strwidth(name[1], cex=text_cex)))
        for(i in seq(along=y)[-1]) {
            for(j in 1:maxy) {
                if(start[i] > maxx[j] + xpad) {
                    y[i] <- j
                    maxx[j] <- max(c(end[i], start[i] + strwidth(name[i], cex=text_cex)))
                    break
                }
            }
            if(is.na(y[i])) { # need new row
                y[i] <- maxy + 1
                maxy <- maxy + 1
                maxx[maxy] <- max(c(end[i], start[i] + strwidth(name[i], cex=text_cex)))
            }
        }

        maxy <- max(c(maxy, minrow))
        height <- 1/maxy

    }

    ypos <- seq(0, by=height, length=maxy)
    y <- ypos[y]

    colors <- rep(colors, length(y))
    for(i in seq(along=start)) {
        rect(start[i], y[i]+(height*padding/4), end[i], y[i]+height/2-(height*padding/4),
             col=colors[i], border=colors[i],
             lend=1, ljoin=1)
        text(start[i], y[i]+height*0.75,
             name[i], adj=c(0, 0.5), col=colors[i],
             cex=text_cex)
    }
}
