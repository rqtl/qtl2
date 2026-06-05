#' Heat map of scan1 results with multiple traits
#'
#' Heat map of multiple scan1 results
#'
#' @param x An object of class `"scan1"`, as output by [scan1()].
#'
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#'
#' @param y Vector with length equal to `ncol(x)`, defining the
#'     y-axis scale of the heat map. (If the multiple LOD scores are
#'     for a quantitative variable like time, these are the
#'     corresponding times.)
#'
#' @param chr Selected chromosomes to plot; a vector of character
#' strings.
#'
#' @param gap Gap between chromosomes. The default is 1% of the total genome length.
#'
#' @param zlim Z-axis limits (a pair of numbers)
#'
#' @param color_scheme Color scheme for the heatmap (ignored if `col` is provided).
#'
#' @param col Optional vector of colors for the heatmap.
#'
#' @param n_colors Number of z-axis colors (ignored if `col` is provided).
#'
#' @param swap_axes If TRUE, put chromosome position on y-axis and lod columns on x-axis
#' (default is to have chromosome position on x-axis and lod columns on y-axis)
#'
#' @param chrlines Color of lines at chromosome breaks (NULL for default choices; NA to skip them).
#'
#' @param rescale_by Rescale the LOD scores by column, chromosome, or both, by subtracting the
#' maximum value. The default is not to rescale.
#'
#' @param ... Additional graphics paramaters.
#'
#' @return None.
#'
#' @seealso [plot_scan1()], [plot_lodpeaks()], [graphics::image()]; [plot_colorscale()]
#'
#' @importFrom graphics image par axis title
#' @importFrom grDevices rainbow heat.colors terrain.colors topo.colors gray
#'
#' @export
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
#' out <- scan1(probs, grav2$pheno)
#' plot_scan1_heatmap(out, map)
#'
#' # make the y-axis quantitative
#' y <- grav2$phenocovar$time
#' plot_scan1_heatmap(out, map, y=y, ylab=colnames(grav2$phenocovar))
#'
#' # plot with colors white -> darkslateblue, and with color scale
#' col <- colorRampPalette(c("white", "darkslateblue"))(256)
#' zlim <- c(0, maxlod(out))
#' layout(cbind(1,2), width=c(3,1))
#' plot_scan1_heatmap(out, map, y=y, ylab=colnames(grav2$phenocovar),
#'                    col=col, zlim=zlim)
#' plot_colorscale(col=col, zlim=zlim)

plot_scan1_heatmap <-
    function(x, map, chr=NULL, y=NULL, gap=NULL, zlim=NULL,
             color_scheme=c("viridis", "gray", "revgray", "heat", "terrain", "topo", "rainbow"),
             col=NULL, n_colors=256, swap_axes=FALSE, chrlines=NULL,
             rescale_by=c("none", "col", "chr", "both"), ...)
{
    if(is.null(map)) stop("map is NULL")

    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    if(!is.matrix(x) && !is.data.frame(x)) stop("x must be a matrix or data frame")

    # if input x is a vector, try turning it into a matrix, since we'll be pull out row names
    if(!is.matrix(x) && !is.data.frame(x) && is.numeric(x)) x <- as.matrix(x)

    # subset chromosomes
    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        x <- subset_scan1(x, map, chr)
        map <- map[chri]
    }

    # align scan1 output and map
    tmp <- align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map

    if(is.null(gap)) gap <- sum(chr_lengths(map))/100
    if(!is_nonneg_number(gap)) stop("gap should be a single non-negative number")

    if(is.null(col)) {
        color_scheme <- match.arg(color_scheme)

        i <- seq(0, 1, length=n_colors)

        if(color_scheme=="viridis") {
            col <- viridis_qtl2(n_colors)
            if(is.null(chrlines)) chrlines <- "white"
        } else if(color_scheme=="gray") {
            col <- grDevices::gray(i)
            if(is.null(chrlines)) chrlines <- "violetred"
        } else if(color_scheme=="revgray") {
            col <- rev(grDevices::gray(i))
            if(is.null(chrlines)) chrlines <- "violetred"
        } else if(color_scheme=="heat") {
            col <- grDevices::heat.colors(n_colors)
            if(is.null(chrlines)) chrlines <- "white"
        } else if(color_scheme=="terrain") {
            col <- grDevices::terrain.colors(n_colors)
            if(is.null(chrlines)) chrlines <- "black"
        } else if(color_scheme=="topo") {
            col <- grDevices::topo.colors(n_colors)
            if(is.null(chrlines)) chrlines <- "white"
        } else if(color_scheme=="rainbow") {
            col <- rev(grDevices::rainbow(n_colors, start=0, end=2/3))
            if(is.null(chrlines)) chrlines <- "white"
        } else stop('color_scheme "', color_scheme, '" not allowed')
    } else {
        if(is.null(chrlines)) chrlines <- "white"
    }

    # rescale the lod scores if desired
    x <- rescale_lod(x, map, rescale_by)

    plot_scan1_heatmap_internal <-
        function(x, map, y=NULL, gap, col, chrlines, swap_axes=FALSE, xlab=NULL, ylab=NULL,
                 xlim=NULL, ylim=NULL, zlim=NULL,
                 mgp=NULL, mgp.x=NULL, mgp.y=NULL,
                 las=1, xaxt="s", yaxt="s", chrlines_lwd=2, ...)
    {
        # get x-axis range
        if(is.null(xlim)) {
            start <- xpos_scan1(map, chr=names(map), gap=gap,
                                names(map)[1], min(map[[1]]))
            end <- xpos_scan1(map, chr=names(map), gap=gap,
                              names(map)[length(map)], max(map[[length(map)]]))

            xlim <- c(start, end)
        }
        if(is.null(zlim)) {
            zlim <- range(as.numeric(x), na.rm=TRUE)
            if(zlim[1] > 0) zlim[1] <- 0
        }

        if(is.null(mgp.x)) {
            if(!is.null(mgp)) mgp.x <- mgp
            else mgp.x <- c(1.7, 0.3, 0)
        }
        if(is.null(mgp.y)) {
            if(!is.null(mgp)) mgp.y <- mgp
            else mgp.y <- c(2.0, 0.4, 0)
        }

        xpos <- unlist(lapply(seq_along(map), function(chr) xpos_scan1(map, chr=names(map), gap=gap, names(map)[chr], map[[chr]])))

        if(any(diff(xpos)==0)) { # deal with identical positions
            xpos <- xpos + seq(0, 1e-6, length=length(xpos))
        }

        if(is.null(y)) {
            y <- 1:ncol(x)
            ylim <- c(0.5, ncol(x)+0.5)
            yaxis_quant <- FALSE
        } else {
            stopifnot(length(y) == ncol(x))
            ylim <- range(y)
            yaxis_quant <- TRUE
        }

        if(swap_axes) {
            graphics::image(y, xpos, t(x), xaxs="i", yaxs="i", xlim=ylim, ylim=xlim, zlim=zlim, xlab="", ylab="", xaxt="n", yaxt="n", col=col)
        } else {
            graphics::image(xpos, y, x, xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, zlim=zlim, xlab="", ylab="", xaxt="n", yaxt="n", col=col)
        }

        # chromosome axis labels
        chr <- names(map)
        if(length(chr) > 1) {
            chr_midpt <- sapply(seq_along(map), function(i) xpos_scan1(map, names(map), gap=gap,
                                                                       names(map)[i], mean(range(map[[i]], na.rm=TRUE))))

            for(i in seq_along(chr_midpt)) {
                if(!swap_axes && xaxt != "n") graphics::axis(side=1, at=chr_midpt[i], chr[i], mgp=mgp.x, tick=FALSE, las=las)
                if(swap_axes && yaxt != "n") graphics::axis(side=2, at=chr_midpt[i], chr[i], mgp=mgp.y, tick=FALSE, las=las)
            }
        } else { # single chromosome; put position on axis
            if(!swap_axes && xaxt != "n") graphics::axis(side=1, at=pretty(xpos, n=7), mgp=mgp.x, tick=FALSE, las=las)
            if(swap_axes && yaxt != "n") graphics::axis(side=2, at=pretty(xpos, n=7), mgp=mgp.y, tick=FALSE, las=las)
        }

        # lod column axis labels
        if(!swap_axes && yaxt != "n") {
            if(yaxis_quant) graphics::axis(side=2, at=pretty(y), mgp=mgp.y, tick=FALSE, las=las)
            else graphics::axis(side=2, at=1:ncol(x), colnames(x), mgp=mgp.y, tick=FALSE, las=las)
        }
        if(swap_axes && xaxt != "n") {
            if(yaxis_quant) graphics::axis(side=1, at=pretty(y), mgp=mgp.x, tick=FALSE, las=las)
            else graphics::axis(side=1, at=1:ncol(x), colnames(x), mgp=mgp.x, tick=FALSE, las=las)
        }

        # axis titles
        if(!is.null(xlab) && xlab != "") graphics::title(xlab=xlab, mgp=mgp.x)
        if(!is.null(ylab) && ylab != "") graphics::title(ylab=ylab, mgp=mgp.y)

        if(!is.na(chrlines) && !is.null(chrlines) && length(map)>1) {
            minpos <- sapply(map, min)
            chrlines_pos <- sapply(seq_along(map)[-1], function(chri) xpos_scan1(map, names(map), gap=gap, names(map)[chri], minpos[chri])-gap/2)
            if(swap_axes) abline(h=chrlines_pos, col=chrlines, lwd=chrlines_lwd)
            else abline(v=chrlines_pos, col=chrlines, lwd=chrlines_lwd)
        }

    }

    plot_scan1_heatmap_internal(x, map, y=y, gap=gap, col=col, chrlines=chrlines, swap_axes=swap_axes, zlim=zlim, ...)

}



# rescale a matrix of lod scores
#
# none: leave as is
# col: subtract off the max in each column
# chr: on each chr, subtract off the overall max
# both: on each chr, subtract off the max in each column
rescale_lod <-
    function(x, map, rescale_by=c("none", "col", "chr", "both"))
{
    rescale_by <- match.arg(rescale_by)

    if(rescale_by=="col") {
        # subtract off the maximum in each column
        x <- t(t(x) - apply(x, 2, max, na.rm=TRUE))
    } else if(rescale_by=="chr") {
        # on each chromosome, subtract off the overall maximum
        for(chr in names(map)) {
            mar <- names(map[[chr]])
            x[mar,] <- x[mar,] - max(x[mar,], na.rm=TRUE)
        }
    } else if(rescale_by=="both") {
        # on each chromosome, subtract off the max in each column
        for(chr in names(map)) {
            mar <- names(map[[chr]])
            z <- x[mar,,drop=FALSE]
            z <- t(t(z) - apply(z, 2, max, na.rm=TRUE))
            x[mar,] <- z
        }
    }

    x
}
