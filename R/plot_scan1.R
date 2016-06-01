#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param x Output of \code{\link[qtl2scan]{scan1}}.
#'
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). Only one value allowed.
#'
#' @param chr Selected chromosomes to plot; a vector of character
#' strings.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param gap Gap between chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @seealso \code{\link{plot_coef}}, \code{\link{plot_coefCC}}, \code{\link{plot_snpasso}}
#'
#' @export
#' @importFrom graphics plot rect lines par axis title abline box
#'
#' @return None.
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
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' library(qtl2scan)
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # plot the results for selected chromosomes
#' ylim <- c(0, maxlod(out)*1.02) # need to strip class to get overall max LOD
#' chr <- c(2,7,8,9,15,16)
#' plot(out, chr=chr, ylim=ylim)
#' plot(out, lodcolumn=2, chr=chr, col="violetred", add=TRUE)
#' legend("topleft", lwd=2, col=c("darkslateblue", "violetred"), colnames(out$lod),
#'        bg="gray90")
#'
#' # plot just one chromosome
#' plot(out, chr=8, ylim=ylim)
#' plot(out, chr=8, lodcolumn=2, col="violetred", add=TRUE)
#'
#' # lodcolumn can also be a column name
#' plot(out, lodcolumn="liver", ylim=ylim)
#' plot(out, lodcolumn="spleen", col="violetred", add=TRUE)
plot_scan1 <-
    function(x, lodcolumn=1, chr=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    # pull out map
    map <- x$map
    if(is.null(map)) stop("No map found in the input")
    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    # pull out lod scores
    if(length(lodcolumn) > 1) { # If length > 1, take first value
        warning("lodcolumn should have length 1; one first element used.")
        lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(x$lod))
        if(is.na(tmp))
            stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(x$lod))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x$lod), ")")
    lod <- x$lod[,lodcolumn]

    # subset chromosomes
    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        map <- map[chri]
        lod <- lod[unlist(lapply(map, names))]
    }

    # internal function; trick to be able to pull things out of "..."
    #    but still have some defaults for them
    plot_scan1_internal <-
        function(map, lod, add=FALSE, gap,
                 bgcolor, altbgcolor,
                 lwd=2, col="darkslateblue", xlab=NULL, ylab="LOD score",
                 xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
                 main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
                 mgp=NULL, las=1,
                 hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
                 vlines=NULL, vlines.col="white", vlines.lwd=1, vlines.lty=1,
                 ...)
        {
            dots <- list(...)
            onechr <- (length(map)==1) # single chromosome

            xpos <- map_to_xpos(map, gap)
            chrbound <- map_to_boundaries(map, gap)

            if(!add) { # new plot
                if(is.null(ylim))
                    ylim <- c(0, max(lod, na.rm=TRUE)*1.02)

                if(is.null(xlim)) {
                    xlim <- range(xpos, na.rm=TRUE)
                    if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
                }

                if(is.null(xlab)) {
                    if(onechr) {
                        if(names(map) == " ") xlab <- "Position"
                        else xlab <- paste("Chr", names(map), "position")
                    }
                    else xlab <- "Chromosome"
                }

                # margin parameters
                if(!is.null(mgp)) mgp.x <- mgp.y <- mgp

                # make basic plot
                plot(xpos, lod, xlab="", ylab="", xlim=xlim, ylim=ylim,
                     xaxs=xaxs, yaxs=yaxs, xaxt="n", yaxt="n", type="n",
                     main=main)

                # add background rectangles
                u <- par("usr")
                if(!is.null(bgcolor))
                    rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
                if(!is.null(altbgcolor) && !onechr) {
                    for(i in seq(2, ncol(chrbound), by=2))
                        rect(chrbound[1,i], u[3], chrbound[2,i], u[4], col=altbgcolor, border=NA)
                }

                # include axis labels?
                if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
                if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

                # add x axis unless par(xaxt="n")
                if(dots$xaxt != "n") {
                    if(onechr) {
                        axis(side=1, at=pretty(xlim), mgp=mgp.x, las=las, tick=FALSE)
                        if(!(length(vlines)==1 && is.na(vlines))) { # if vlines==NA, skip lines
                            if(is.null(vlines)) vlines <- pretty(xlim)
                            abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
                        }
                    }
                    else {
                        loc <- colMeans(chrbound)
                        odd <- seq(1, length(map), by=2)
                        even <- seq(2, length(map), by=2)
                        axis(side=1, at=loc[odd], names(map)[odd],
                             mgp=mgp.x, las=las, tick=FALSE)
                        axis(side=1, at=loc[even], names(map)[even],
                             mgp=mgp.x, las=las, tick=FALSE)
                    }
                }
                else if(onechr) { # if dots$xaxt="n" but one chr and vlines explicit, plot them
                    if(!is.null(vlines) && !(length(vlines)==1 && is.na(vlines)))
                        abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
                }

                # add y axis unless par(yaxt="n")
                if(dots$yaxt != "n") {
                    axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
                    if(!(length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
                        if(is.null(hlines)) hlines <- pretty(ylim)
                        abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
                    }
                }
                dots$xaxt <- dots$yaxt <- NULL # delete those

                # x and y axis labels
                title(xlab=xlab, mgp=mgp.x)
                title(ylab=ylab, mgp=mgp.y)
            }

            # plot each chromosome
            indexes <- map_to_index(map)
            for(i in seq(along=indexes))
                lines(xpos[indexes[[i]]], lod[indexes[[i]]],
                           lwd=lwd, col=col, ...)

            # add box just in case
            box()
        }

    # make the plot
    plot_scan1_internal(map=map, lod=lod, add=add, gap=gap,
                       bgcolor=bgcolor, altbgcolor=altbgcolor,
                       ...)
}


#' @export
#' @rdname plot_scan1
plot.scan1 <-
    function(x, lodcolumn=1, chr=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    # if snp asso result, use plot_snpasso() with just reduced snps; otherwise defaults
    if(!is.null(x$snpinfo)) {
        plot_snpasso(x, add=add, gap=gap, bgcolor=bgcolor,
                     altbgcolor=altbgcolor, ...)
    }
    else { # mostly, use plot_scan1()
        plot_scan1(x, lodcolumn=lodcolumn, chr=chr, add=add, gap=gap,
                   bgcolor=bgcolor, altbgcolor=altbgcolor, ...)
    }
}


# convert map to x-axis positions for plot_scan1
map_to_xpos <-
    function(map, gap)
{
    if(length(map)==1) return(map[[1]])

    chr_range <- vapply(map, range, c(0,1), na.rm=TRUE)

    result <- map[[1]]-chr_range[1,1] + gap/2
    for(i in 2:length(map)) {
        result <- c(result,
                    map[[i]] - chr_range[1,i] + gap + max(result, na.rm=TRUE))
    }
    result
}

# boundaries of chromosomes in plot_scan1
# first row: left edges
# second row: right edges
map_to_boundaries <-
    function(map, gap)
{
    if(length(map)==1)
        return(cbind(range(map[[1]], na.rm=TRUE)))

    # range of each chromosome
    chr_range <- lapply(map, range, na.rm=TRUE)

    # corresponding xpos, as matrix with two rows
    startend <- matrix(map_to_xpos(chr_range, gap), nrow=2)

    startend[1,] <- startend[1,] - gap/2
    startend[2,] <- startend[2,] + gap/2

    startend
}

# convert map to list of indexes to LOD vector
map_to_index <-
    function(map)
{
    if(length(map)==1) {
        map[[1]] <- seq(along=map[[1]])
        return(map)
    }

    lengths <- vapply(map, length, 0)
    split(1:sum(lengths), rep(seq(along=map), lengths))
}
