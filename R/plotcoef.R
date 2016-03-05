#' Plot QTL effects along chromosome
#'
#' Plot estimated QTL effects along a chromosomes.
#'
#' @param coef Estimated QTL effects ("coefficients") as obtained from
#' \code{\link[qtl2scan]{scan1coef}} or
#' \code{\link[qtl2scan]{scan1coef_lmm}}.
#'
#' @param columns Vector of columns to plot
#'
#' @param col Vector of colors, same length as \code{columns}.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param gap Gap between chromosomes.
#'
#' @param ylim y-axis limits
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @export
#'
#' @details
#' \code{plotcoefCC()} is the same as \code{plotcoef()}, but forcing
#' \code{columns=1:8} and using the Collaborative Cross colors,
#' \code{\link{CCcolors}}.
#'
#' @seealso \code{\link{CCcolors}}
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' library(qtl2scan)
#' coef <- scan1coef(probs[,7], pheno, covar)
#'
#' # plot QTL effects
#' plotcoef(coef, columns=1:3, col=c("slateblue", "violetred", "green3"))
plotcoef <-
    function(coef, columns, col, add=FALSE, gap=25, ylim,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    if(missing(columns) || is.null(columns))
        columns <- 1:ncol(coef)

    if(missing(col) || is.null(col)) {
        n_col <- length(columns)
        if(n_col == 1) col <- c("slateblue")
        else if(n_col == 2) col <- c("slateblue", "violetred")
        else if (n_col==3) col=c("slateblue", "violetred", "green3")
        else if(n_col <= 8) col <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                "#66A61E", "#E6AB02", "#A6761D", "#666666")[1:n_col]
        else
            stop("With >8 coefficients, you need to provide colors via col")
    } else {
        if(length(col) < length(columns))
            stop("Need at least ", length(columns), " colors")
    }

    map <- attr(coef, "map")
    if(is.null(map)) stop("Input needs to contain a map attribute")

    if(missing(ylim) || is.null(ylim)) {
        ylim <- range(coef[,columns], na.rm=TRUE)
        d <- diff(ylim) * 0.02 # add 2% on either side
        ylim <- ylim + c(-d, d)
    }

    plotscan1(coef, lodcolumn=columns[1], ylim=ylim, col=col[1], add=add,
              gap=gap, bgcolor=bgcolor, altbgcolor=altbgcolor, ...)
    if(length(columns) > 1) {
        for(i in seq(along=columns)[-1])
            plotscan1(coef, lodcolumn=columns[i], col=col[i], gap=gap,
                      add=TRUE, ...)
    }
}

#' @export
#' @rdname plotcoef
plotcoefCC <-
    function(coef, add=FALSE, gap=25, ylim=NULL,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    plotcoef(coef, columns=1:8, col=qtl2plot::CCcolors, add=add, gap=gap,
             ylim=ylim, bgcolor=bgcolor, altbgcolor=altbgcolor, ...)
}
