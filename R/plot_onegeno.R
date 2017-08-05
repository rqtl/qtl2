#' Plot one individual's genome-wide genotypes
#'
#' Plot one individual's genome-wide genotypes
#'
#' @param geno Imputed phase-known genotypes, as a list of matrices
#'     (as produced by \code{\link[qtl2geno]{maxmarg}}) or a list of
#'     three-dimensional arrays (as produced by \code{\link[qtl2geno]{guess_phase}}).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosomes to plot; a vector of character strings.
#' @param col Vector of colors for the different genotypes.
#' @param na_col Color for missing segments.
#' @param border Color of outer border around chromosome rectangles.
#' @param shift If TRUE, shift the chromosomes so they all start at 0.
#' @param bgcolor Color for background rectangle
#' @param chrwidth Total width of rectangles for each chromosome, as a
#'     fraction of the distance between them.
#' @param ... Additional graphics parameters
#'
#' @export
#' @importFrom graphics plot rect par axis title abline box
plot_onegeno <-
    function(geno, map, ind=1, chr=NULL, col=NULL, na_col="white",
             border="black", shift=FALSE, bgcolor="gray90",
             chrwidth=0.5, ...)
{
    if(is.null(map)) stop("map is NULL")

    stopifnot(chrwidth > 0 && chrwidth < 1)

    # ignore class of geno object
    geno <- unclass(geno)

    # drop all but the target individual
    if(length(ind)>1) {
        ind <- ind[1]
        warning("Only using the first individual")
    }
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

    # shift map to start at 0
    if(shift) map <- lapply(map, function(a) a-min(a,na.rm=TRUE))

    plot_onegeno_internal <-
        function(geno, map, col=NULL, na_col="white",
                 border="black", bgcolor="gray90",
                 chrwidth=0.5,
                 xlab="Chromosome", ylab="Position (Mbp)",
                 ylim=NULL, main="", las=1, xaxs="i", yaxs="r",
                 mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
                 hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
                 vlines=NULL, vlines.col="gray80", vlines.lwd=1, vlines.lty=1,
                 ...)
    {
        dots <- list(...)
        if(is.null(ylim)) ylim <- rev(range(unlist(map), na.rm=TRUE))

        nchr <- length(map)
        xlim <- c(0.5, nchr+0.5)

        # margin parameters
        if(!is.null(mgp)) mgp.x <- mgp.y <- mgp

        plot(0, 0, type="n", xlab="", ylab="", main=main,
             xaxs=xaxs, yaxs=yaxs,
             xaxt="n", yaxt="n", xlim=xlim, ylim=ylim)
        u <- par("usr")
        if(!is.null(bgcolor))
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)

        # include axis labels?
        if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
        if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

        # add x axis unless par(xaxt="n")
        if(dots$xaxt != "n") {
            odd <- seq(1, nchr, by=2)
            axis(side=1, at=odd, names(map)[odd],
                 mgp=mgp.x, las=las, tick=FALSE)
            if(nchr > 1) {
                even <- seq(2, nchr, by=2)
                axis(side=1, at=even, names(map)[even],
                     mgp=mgp.x, las=las, tick=FALSE)
            }
        }
        # add y axis unless par(yaxt="n")
        if(dots$yaxt != "n") {
            axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
        }

        # grid lines
        if(!(length(vlines)==1 && is.na(vlines))) {
            if(is.null(vlines)) vlines <- 1:nchr
            abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
        }
        if(!(length(hlines)==1 && is.na(hlines))) {
            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
        }

        # x and y axis labels
        title(xlab=xlab, mgp=mgp.x)
        title(ylab=ylab, mgp=mgp.y)

        max_geno <- max(unlist(geno), na.rm=TRUE)
        if(is.null(col)) {
            data(CCcolors)
            if(max_geno <= 8) {
                col <- CCcolors
            }
            else {
                warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
                col <- rep(CCcolors, max_geno)
            }
        }
        else if(max_geno > length(col)) {
            warning("not enough colors; recycling them")
            col <- rep(col, max_geno)
        }

        for(i in 1:nchr) {
            g <- geno[[i]]

            # if completely missing the second chr, treat as if we have just the one
            this_chrwidth <- chrwidth
            if(!is.matrix(g) && all(is.na(g[,,2]))) {
                g <- rbind(g[,,1]) # make it a row matrix
                this_chrwidth <- this_chrwidth/2
            }

            if(is.matrix(g)) { # phase-known
                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i+this_chrwidth/2, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)

                addgenorect(g[1,], map[[i]], i-this_chrwidth/2, i+this_chrwidth/2,
                            col=col)

                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i+this_chrwidth/2, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
            }
            else {
                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)

                addgenorect(g[1,,1], map[[i]], i-this_chrwidth/2, i,
                            col=col)
                addgenorect(g[1,,2], map[[i]], i+this_chrwidth/2, i,
                            col=col)

                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
            }
        }

        box()
    }

    plot_onegeno_internal(geno, map, col=col, na_col=na_col, border=border,
                          bgcolor=bgcolor, chrwidth=chrwidth, ...)


}

# add rectangles for the genotypes
addgenorect <-
    function(geno, map, x1, x2, col)
{
    intervals <- geno2intervals(geno, map)
    if(is.null(intervals) || nrow(intervals) < 1) return(NULL)

    for(i in 1:nrow(intervals)) {
        rect(x1, intervals[i,1],
             x2, intervals[i,2],
             col=col[intervals[i,3]],
             border=NA, lend=1, ljoin=1)
    }
}


# convert vector of integer genotypes to intervals with common genotypes
# (start, end, genotype)
geno2intervals <-
    function(geno, map)
{
    if(all(is.na(geno))) return(NULL)

    stopifnot(length(geno) == length(map))

    # drop missing values
    map <- map[!is.na(geno)]
    geno <- geno[!is.na(geno)]

    d <- diff(geno)
    xo_int <- which(d != 0)

    data.frame(lo=map[c(1,xo_int+1)],
               hi=map[c(xo_int, length(map))],
               geno=geno[c(xo_int, length(map))])

}
