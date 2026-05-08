#' Heat map color scale
#'
#' Plot heat map color scale
#'
#' @param x Values that define the color limits (matrix, data frame, or vector); ignored if `zlim` is provided.
#'
#' @param zlim Z-axis limits (a pair of numbers)
#'
#' @param color_scheme Color scheme for the heatmap (ignored if `col` is provided).
#'
#' @param col Optional vector of colors for the heatmap.
#'
#' @param n_colors Number of z-axis colors (ignored if `col` is provided).
#'
#' @param gridlines Color of gridlines (NULL for default choices; NA to skip them).
#'
#' @param swap_axes If TRUE, make the scale horizontal
#'
#' @param ... Additional graphics paramaters.
#'
#' @return None.
#'
#' @seealso [plot_scan1_heatmap()]
#'
#' @importFrom graphics image par axis
#' @importFrom grDevices rainbow heat.colors terrain.colors topo.colors gray
#'
#' @export
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
#' out <- scan1(probs, grav2$pheno)
#'
#' layout(cbind(1,2), widths=c(4, 1))
#' par(mar=c(4.1, 4.1, 1.1, 1.1))
#' plot_scan1_heatmap(out, map, xlab="Chromosome", ylab="Phenotype", mgp.y=c(2.6, 0.3, 0))
#' par(mar=c(15.1, 3.1, 15.1, 1.1))
#' plot_colorscale(out)

plot_colorscale <-
    function(x=NULL, zlim=NULL,
             color_scheme=c("viridis", "gray", "revgray", "heat", "terrain", "topo", "rainbow"),
             col=NULL, n_colors=256, gridlines=NULL, swap_axes=FALSE, ...)
{
    if(is.null(zlim)) {
        if(!is.null(x)) {
            x <- as.numeric(unlist(x))
            zlim <- c(0, max(x, na.rm=TRUE))
        } else {
            stop("Provide either x or zlim")
        }
    }

    if(is.null(col)) {
        color_scheme <- match.arg(color_scheme)

        i <- seq(0, 1, length=n_colors)

        if(color_scheme=="viridis") {
            col <- viridis_qtl2(n_colors)
            if(is.null(gridlines)) gridlines <- "white"
        } else if(color_scheme=="gray") {
            col <- grDevices::gray(i)
            if(is.null(gridlines)) gridlines <- "violetred"
        } else if(color_scheme=="revgray") {
            col <- rev(grDevices::gray(i))
            if(is.null(gridlines)) gridlines <- "violetred"
        } else if(color_scheme=="heat") {
            col <- grDevices::heat.colors(n_colors)
            if(is.null(gridlines)) gridlines <- "white"
        } else if(color_scheme=="terrain") {
            col <- grDevices::terrain.colors(n_colors)
            if(is.null(gridlines)) gridlines <- "black"
        } else if(color_scheme=="topo") {
            col <- grDevices::topo.colors(n_colors)
            if(is.null(gridlines)) gridlines <- "white"
        } else if(color_scheme=="rainbow") {
            col <- rev(grDevices::rainbow(n_colors, start=0, end=2/3))
            if(is.null(gridlines)) gridlines <- "white"
        } else stop('color_scheme "', color_scheme, '" not allowed')
    } else {
        if(is.null(gridlines)) gridlines <- "white"
    }

    plot_colorscale_internal <-
        function(zlim, col, swap_axes=FALSE,
                 gridlines=NULL, gridlines_lwd=1,
                 mgp=NULL, mgp.x=NULL, mgp.y=NULL,
                 las=1, xaxt="s", yaxt="s",
                 xlab="", ylab="", ...)
    {
        xpos <- 1
        ypos <- seq(zlim[1], zlim[2], length=length(col))
        z <- as.matrix(ypos)

        if(swap_axes) {
            graphics::image(ypos, xpos, z, xaxs="i", yaxs="i", zlim=zlim, xlab="", ylab="", xaxt="n", yaxt="n", col=col)
        } else {
            graphics::image(xpos, ypos, t(z), xaxs="i", yaxs="i", zlim=zlim, xlab="", ylab="", xaxt="n", yaxt="n", col=col)
        }

        if(is.null(mgp.x)) {
            if(!is.null(mgp)) mgp.x <- mgp
            else mgp.x <- c(1.7, 0.3, 0)
        }
        if(is.null(mgp.y)) {
            if(!is.null(mgp)) mgp.y <- mgp
            else mgp.y <- c(2.0, 0.4, 0)
        }

        # lod column axis labels
        if(!swap_axes && yaxt != "n") {
            graphics::axis(side=2, at=pretty(ypos), mgp=mgp.y, tick=FALSE, las=las)
            abline(h=pretty(ypos), col=gridlines, lwd=gridlines_lwd)
        }
        if(swap_axes && xaxt != "n") {
            graphics::axis(side=1, at=pretty(ypos), mgp=mgp.x, tick=FALSE, las=las)
            abline(v=pretty(ypos), col=gridlines, lwd=gridlines_lwd)
        }

        # axis titles
        if(!is.null(xlab) && xlab != "") graphics::title(xlab=xlab, mgp=mgp.x)
        if(!is.null(ylab) && ylab != "") graphics::title(ylab=ylab, mgp=mgp.y)
    }

    plot_colorscale_internal(zlim=zlim, col=col, swap_axes=swap_axes, gridlines=gridlines, ...)

}
