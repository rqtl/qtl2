#' Plot multiple individuals' genome-wide genotypes
#'
#' Plot multiple individuals' genome-wide genotypes
#'
#' @param geno Imputed phase-known genotypes, as a list of matrices
#'     (as produced by [maxmarg()]) or a list of
#'     three-dimensional arrays (as produced by [guess_phase()]).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individuals to plot, either a numeric indexes or IDs.
#' @param chr Selected chromosomes to plot; a vector of character strings.
#' @param gap Gap between chromosomes
#' @param col Vector of colors for the different genotypes.
#' @param na_col Color for missing segments.
#' @param chrlines Color for lines separating chromosomes
#' @param swap_axes If TRUE, swap the axes, so that the chromosomes run horizontally.
#' @param ... Additional graphics parameters
#'
#' @return None.
#'
#' @seealso [plot_onegeno()], [plot_genoprob()]
#'
#' @examples
#' # load data and calculate genotype probabilities
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' iron <- iron[1:50, ] # subset to first 50 individuals
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # infer genotypes, as those with maximal marginal probability
#' m <- maxmarg(pr, minprob=0.5)
#'
#' # re-code the X chr, (5,6) -> (1,3)
#' m[["X"]] <- (m[["X"]] - 5)*2 + 1
#'
#' # plot phased genotypes
#' plot_geno(m, map, col=c("#FFDC00", "#00C800", "#0064C9"))
#'
#' # this is more interesting for Diversity Outbred mouse data
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' # subset to first 25 individuals
#' DOex <- DOex[1:25, ]
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#'
#' # infer genotypes, as those with maximal marginal probability
#' m <- maxmarg(pr, minprob=0.5)
#' # guess phase
#' ph <- guess_phase(DOex, m)
#'
#' # plot phased genotypes
#' plot_geno(ph, DOex$gmap)
#' }
#'
#' @export
#' @importFrom graphics image
plot_geno <-
    function(geno, map, ind=NULL, chr=NULL, gap=0,
             col=NULL, na_col="white", chrlines="white",
             swap_axes=FALSE, ...)

{
    if(is.null(geno)) stop("geno is NULL")
    if(is.null(map)) stop("map is NULL")

    # ignore class of geno object
    geno <- unclass(geno)

    # subset individuals
    if(!is.null(ind)) {
        for(i in seq_along(geno)) {
            if(is.matrix(geno[[i]]))
                geno[[i]] <- geno[[i]][ind,,drop=FALSE]
            else {
                if(!is.array(geno[[i]]) || length(dim(geno[[i]])) != 3 ||
                   dim(geno[[i]])[3] != 2)
                    stop("geno should be an array individuals x positions x 2 haplotypes")
                geno[[i]] <- geno[[i]][ind,,,drop=FALSE]
            }
        }
    }

    # find common chr
    chr_geno <- names(geno)
    chr_map <- names(map)
    common_chr <- chr_geno[chr_geno %in% chr_map]
    if(length(common_chr) == 0)
        stop("No chr in common between geno and map")
    geno <- geno[common_chr]
    map <- map[common_chr]

    # subset chr if necessary
    if(!is.null(chr)) {
        if(any(!(chr %in% common_chr))) {
            chr <- chr[chr %in% common_chr]
            if(length(chr) == 0)
                stop("No chromosomes in common between geno, map, and chr")
            warning("Dropping some chr not found in geno and/or map")
        }
        geno <- geno[chr]
        map <- map[chr]
    }

    # same numbers of markers?
    nmar_map <- sapply(map, length)
    nmar_geno <- sapply(geno, ncol)
    if(any(nmar_geno != nmar_map))
        stop("Mismatch between numbers of markers between geno and map on chr ",
             paste(names(geno)[nmar_geno != nmar_map], collapse=", "))
    for(i in seq_along(geno)) {
        if(!all(names(map[[i]]) == colnames(geno[[i]])))
            stop("Mismatch between marker names on chr ", names(geno)[i])
    }

    if(length(dim(geno[[1]]))>2) { # genotypes are 2d arrays; stack them
        geno <- lapply(geno, function(a) {
            d3 <- dim(a)[3]
            result <- matrix(nrow=nrow(a)*d3, ncol=ncol(a))
            rownames(result) <- 1:nrow(result)
            colnames(result) <- colnames(geno)
            for(i in 1:d3) {
                result[(1:nrow(a))*d3-(d3-i),] <- a[,,i]
                rownames(result)[(1:nrow(a))*d3-(d3-i)] <- rownames(a)
            }
            result
        })
    }

    plot_geno_internal <-
        function(geno, map, col=NULL, na_col="white",
                 chrlines="white", chrlines_lwd=2,
                 swap_axes=FALSE,
                 xlab=NULL, ylab=NULL,
                 xlim=NULL, ylim=NULL, las=1,
                 mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
                 xaxt="s", yaxt="s",
                 ...)
    {
        # margin parameters
        if(!is.null(mgp)) mgp.x <- mgp.y <- mgp

        # paste genos together
        geno <- do.call("cbind", geno)

        # z-axis range and breaks
        geno[is.na(geno)] <- 0
        n_colors <- max(geno)+1
        zlim <- c(0, max(geno))
        breaks <- seq(-0.5, n_colors-0.5, by=1)

        if(is.null(col)) {
            if(n_colors==3) col <- c("#FFDC00", "#0064C9")
            else if(n_colors==4) col <- c("#FFDC00", "#00C800", "#0064C9")
            else col <- qtl2::CCcolors
        }
        col <- c(na_col, col)
        if(length(breaks) > length(col)+1)
            stop("Need more colors: at least ", n_colors-1)
        col <- col[1:(length(breaks)-1)]

        # get x-axis range
        if(is.null(xlim)) {
            start <- xpos_scan1(map, chr=names(map), gap=gap,
                                names(map)[1], min(map[[1]]))
            end <- xpos_scan1(map, chr=names(map), gap=gap,
                              names(map)[length(map)], max(map[[length(map)]]))

            xlim <- c(start, end)
        }
        if(is.null(ylim)) {
            ylim <- c(0.5, nrow(geno)+0.5)
        }

        xpos <- unlist(lapply(seq_along(map), function(chr) xpos_scan1(map, chr=names(map), gap=gap, names(map)[chr], map[[chr]])))

        if(any(diff(xpos)==0)) { # deal with identical positions
            xpos <- xpos + seq(0, 1e-6, length=length(xpos))
        }

        ypos <- 1:nrow(geno)

        if(!swap_axes) image(xpos, ypos, t(geno), breaks=breaks, col=col, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
        if(swap_axes)  image(ypos, xpos, geno,    breaks=breaks, col=col, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n")

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

        # ind column axis labels
        if(!swap_axes && yaxt != "n") graphics::axis(side=2, at=1:nrow(geno), rownames(geno), mgp=mgp.y, tick=FALSE, las=las)
        if(swap_axes && xaxt != "n") graphics::axis(side=1, at=1:nrow(geno), rownames(geno), mgp=mgp.x, tick=FALSE, las=las)

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

    plot_geno_internal(geno, map, gap=gap, col=col, na_col=na_col,
                       swap_axes=swap_axes, chrlines=chrlines, ...)


}
