#' Plot QTL effects along chromosome
#'
#' Plot estimated QTL effects along a chromosomes.
#'
#' @md
#'
#' @param x Estimated QTL effects ("coefficients") as obtained from
#' [scan1coef()].
#'
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#'
#' @param columns Vector of columns to plot
#'
#' @param col Vector of colors, same length as `columns`. If
#' NULL, some default choices are made.
#'
#' @param scan1_output If provided, we make a two-panel plot with
#' coefficients on top and LOD scores below. Should have just one LOD
#' score column; if multiple, only the first is used.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param gap Gap between chromosomes.
#'
#' @param top_panel_prop If `scan1_output` provided, this gives the
#' proportion of the plot that is devoted to the top panel.
#'
#' @param legend Location of legend, such as `"bottomleft"` or `"topright"` (NULL for no legend)
#'
#' @param ... Additional graphics parameters.
#'
#' @export
#' @importFrom graphics layout par
#'
#' @details
#' `plot_coefCC()` is the same as `plot_coef()`, but forcing
#' `columns=1:8` and using the Collaborative Cross colors,
#' [CCcolors].
#'
#' @seealso [CCcolors], [plot_scan1()], [plot_snpasso()]
#'
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color, and things
#' like `ylab` and `ylim`. These are not included as formal
#' parameters in order to avoid cluttering the function definition.
#'
#' In the case that `scan1_output` is provided, `col`,
#' `ylab`, and `ylim` all control the panel with estimated
#' QTL effects, while `col_lod`, `ylab_lod`, and
#' `ylim_lod` control the LOD curve panel.
#'
#' If `legend` is indicated so that a legend is shown, `legend_lab`
#' controls the labels in the legend, and `legend_ncol` indicates the
#' number of columns in the legend.
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' coef <- scan1coef(probs[,7], pheno, addcovar=covar)
#'
#' # plot QTL effects (note the need to subset the map object, for chromosome 7)
#' plot(coef, map[7], columns=1:3, col=c("slateblue", "violetred", "green3"))
plot_coef <-
    function(x, map, columns=NULL, col=NULL, scan1_output=NULL,
             add=FALSE, gap=25, top_panel_prop=0.65,
             legend=NULL, ...)
{
    if(is.null(map)) stop("map is NULL")

    # align scan1 output and map
    tmp <- align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map

    if(nrow(x) != length(unlist(map)))
        stop("nrow(x) [", nrow(x), "] != number of positions in map [",
             length(unlist(map)), "]")

    if(!is.null(scan1_output)) { # call internal function for both coef and LOD
        return(plot_coef_and_lod(x, map, columns=columns, col=col, scan1_output=scan1_output,
                                 gap=gap, xaxt=NULL, top_panel_prop=top_panel_prop,
                                 legend=legend, ...))
    }

    if(is.null(columns))
        columns <- seq_len(ncol(x))

    if(is.null(col)) {
        n_col <- length(columns)
        if(n_col == 1) col <- c("slateblue")
        else if(n_col == 2) col <- c("slateblue", "violetred")
        else if (n_col==3) col=c("slateblue", "violetred", "green3")
        else if(n_col <= 8) col <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                "#66A61E", "#E6AB02", "#A6761D", "#666666")[seq_len(n_col)]
        else
            stop("With >8 coefficients, you need to provide colors via col")
    } else {
        if(length(col) < length(columns))
            stop("Need at least ", length(columns), " colors")
    }

    plot_coef_internal <-
        function(x, map, columns, ylim=NULL, col, add, gap, bgcolor="gray90", altbgcolor="gray85",
                 ylab="QTL effects", legend_lab=NULL, legend_ncol=NULL, ...)
        {
            if(is.null(ylim)) {
                ylim <- range(unclass(x)[,columns], na.rm=TRUE)
                d <- diff(ylim) * 0.02 # add 2% on either side
                ylim <- ylim + c(-d, d)
            }

            plot_scan1(x, map, lodcolumn=columns[1], ylim=ylim, col=col[1], add=add,
                       gap=gap, bgcolor=bgcolor, altbgcolor=altbgcolor,
                       ylab=ylab, ...)
            if(length(columns) > 1) {
                for(i in seq(along=columns)[-1])
                    plot_scan1(x, map, lodcolumn=columns[i], col=col[i], gap=gap,
                               add=TRUE, ...)
            }

            if(!is.null(legend)) {
                if(is.null(legend_lab)) {
                    if(!is.null(names(col))) legend_lab <- names(col)
                    else legend_lab <- colnames(x)[columns]
                }
                if(length(legend_lab) > length(columns))
                    legend_lab <- legend_lab[seq_along(columns)]
                this_col <- rep(col, length(legend_lab))[seq_along(legend_lab)]
                if(is.null(legend_ncol)) legend_ncol <- ifelse(length(legend_lab) > 4, 2, 1)
                legend(legend, lwd=2, col=this_col, legend_lab, bg=bgcolor,
                       ncol=legend_ncol)
            }

        }
    plot_coef_internal(unclass(x), map, columns=columns, col=col, add=add, gap=gap, ...)
}

#' @export
#' @rdname plot_coef
plot_coefCC <-
    function(x, map, columns=1:8, scan1_output=NULL, add=FALSE, gap=25,
             top_panel_prop=0.65, legend=NULL, ...)
{
    plot_coef(x, map, columns=columns, col=qtl2::CCcolors[columns],
              scan1_output=scan1_output, add=add, gap=gap,
              top_panel_prop=top_panel_prop, legend=legend, ...)
}

#' @export
#' @rdname plot_coef
plot.scan1coef <-
    function(x, map, columns=1, col=NULL, scan1_output=NULL, add=FALSE, gap=25,
             top_panel_prop=0.65, legend=NULL, ...)
{
    plot_coef(x, map, columns=columns, col=col, scan1_output=scan1_output,
              add=add, gap=gap, top_panel_prop=top_panel_prop, legend=legend,
              ...)
}
