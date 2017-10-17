#' Plot phenotype vs genotype
#'
#' Plot phenotype vs genotype for a single putative QTL and a single phenotype.
#'
#' @param geno Vector of genotypes, as produced by
#' \code{\link[qtl2geno]{maxmarg}} with specific \code{chr} and
#' \code{pos}.
#' @param pheno Vector of phenotypes.
#' @param sort If TRUE, sort genotypes from largest to smallest.
#' @param SEmult If specified, interval estimates of the within-group
#' averages will be displayed, as \code{mean +/- SE * SEmult}.
#' @param pooledSD If TRUE and \code{SEmult} is specified, calculated
#' a pooled within-group SD. Otherwise, get separate estimates of
#' the within-group SD for each group.
#' @param jitter Amount to jitter the points horizontally, if a vector
#' of length > 0, it is taken to be the actual jitter amounts
#' (with values between -0.5 and 0.5).
#' @param seg_width Width of segments at the estimated within-group averages
#' @param seg_lwd Line width used to plot estimated within-group averages
#' @param seg_col Line color used to plot estimated within-group averages
#' @param bgcolor Background color for the plot.
#' @param hlines Locations of horizontal grid lines.
#' @param hlines_col Color of horizontal grid lines
#' @param hlines_lty Line type of horizontal grid lines
#' @param hlines_lwd Line width of horizontal grid lines
#' @param vlines_col Color of vertical grid lines
#' @param vlines_lty Line type of vertical grid lines
#' @param vlines_lwd Line width of vertical grid lines
#' @param force_labels If TRUE, force all genotype labels to be shown.
#' @param alternate_labels If TRUE, place genotype labels in two rows
#' @param omit_points If TRUE, omit the points, just plotting the averages (and, potentially, the +/- SE intervals).
#' @param ... Additional graphics parameters, passed to \code{\link[graphics]{plot}}.
#'
#' @export
#' @importFrom graphics par plot segments title axis points
#' @importFrom stats lm runif sd
#'
#' @seealso \code{\link{plot_coef}}
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # inferred genotype at a 28.6 cM on chr 16
#' geno <- maxmarg(probs, map, chr=16, pos=28.6, return_char=TRUE)
#'
#' # plot phenotype vs genotype
#' plot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)))
plot_pxg <-
    function(geno, pheno, sort=TRUE, SEmult=NULL, pooledSD=TRUE,
             jitter=0.2, bgcolor="gray90",
             seg_width=0.4, seg_lwd=2, seg_col="black",
             hlines=NULL, hlines_col="white", hlines_lty=1, hlines_lwd=1,
             vlines_col="gray80", vlines_lty=1, vlines_lwd=3,
             force_labels=TRUE, alternate_labels=FALSE,
             omit_points=FALSE, ...)
{
    if(length(jitter) > 1 && length(jitter) != length(geno)) {
        stop("length(jitter) > 1 but length(jitter) != length(geno)")
        if(any(!is.na(jitter) & (jitter < -0.5 | jitter > 0.5)))
            stop("jitter values should be in [-0.5, +0.5]")
    }
    if(length(jitter) == 1) {
        if(jitter < 0 || jitter > 0.5)
            stop("jitter should be in [0,0.5]")
        jitter <- runif(length(geno), -jitter, jitter)
    }
    names(jitter) <- names(geno)

    # align geno and pheno
    ind_g <- names(geno)
    ind_p <- names(pheno)
    if(!is.null(ind_g) && !is.null(ind_p)) {
        ind <- ind_g[ind_g %in% ind_p]
        if(length(ind) == 0)
            stop("No individuals in common between geno and pheno")
        geno <- geno[ind]
        jitter <- jitter[ind]
        pheno <- pheno[ind]
    }
    else {
        if(length(geno) != length(pheno))
            stop("length(geno) != length(pheno)")
    }

    # genotype groups
    if(!is.factor(geno)) geno <- factor(geno)
    ugeno <- levels(geno)
    me <- tapply(pheno, geno, mean, na.rm=TRUE)
    if(sort) {
        o <- order(me, decreasing=TRUE, na.last=TRUE)
        ugeno <- names(me[o])
        geno <- factor(geno, levels=ugeno)
        me <- tapply(pheno, geno, mean, na.rm=TRUE)
    }

    # calculate intervals
    se <- NULL
    if(!is.null(SEmult)) {
        if(SEmult < 0) stop("SEmult should be >= 0")
        if(pooledSD) {
            sigma <- summary(lm(pheno ~ geno))$sigma
            se <- sigma / sqrt(tapply(pheno, geno, function(a) sum(!is.na(a))))
        }
        else {
            se <- tapply(pheno, geno, function(a) sd(a, na.rm=TRUE)/sqrt(sum(!is.na(a))))
        }
    }

    plot_pxg_internal <-
        function(geno, pheno, bgcolor="gray90",
                 seg_width=0.2, seg_lwd=2, seg_col="slateblue",
                 hlines=NULL, hlines_col="white", hlines_lty=1, hlines_lwd=1,
                 vlines_col="gray80", vlines_lty=1, vlines_lwd=3,
                 xlim=c(0.5, length(ugeno)+0.5), ylim=NULL,
                 xaxs="i", yaxs="r", xlab="Genotype", ylab="Phenotype",
                 mgp=c(2.6, 0.3, 0), mgp.x=mgp, mgp.y=mgp, las=1,
                 pch=21, bg="lightblue", ...)
        {
            # get y-axis limits
            if(is.null(ylim)) {
                if(omit_points) {
                    ylim <- range(me, na.rm=TRUE)
                } else {
                    ylim <- range(pheno, na.rm=TRUE)
                }
                if(!is.null(SEmult)) {
                    ylim <- range(c(ylim, me-se*SEmult, me+se*SEmult), na.rm=TRUE)
                }
            }

            plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n",
                 yaxt="n", xaxs=xaxs, yaxs=yaxs)
            title(xlab=xlab, mgp=mgp.x)
            title(ylab=ylab, mgp=mgp.y)

            u <- par("usr")
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)

            x <- seq_along(ugeno)
            abline(v=x, lwd=vlines_lwd, lty=vlines_lty, col=vlines_col)
            if(force_labels) {
                if(alternate_labels) {
                    for(i in x)
                        axis(side=1, at=i, labels=ugeno[i], mgp=mgp.x, tick=FALSE,
                             line=2 - (i %% 2)-1)
                } else {
                    for(i in x)
                        axis(side=1, at=i, labels=ugeno[i], mgp=mgp.x, tick=FALSE)
                }
            } else {
                if(alternate_labels) {
                    odd <- seq(1, length(ugeno), by=2)
                    axis(side=1, at=x[odd], labels=ugeno[odd], mgp=mgp.x, tick=FALSE)
                    if(length(ugeno) > 1) {
                        axis(side=1, at=x[-odd], labels=ugeno[-odd], mgp=mgp.x, tick=FALSE,
                             line=1)
                    }
                } else {
                    axis(side=1, at=x, labels=ugeno, mgp=mgp.x, tick=FALSE)
                }
            }

            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, lwd=hlines_lwd, lty=hlines_lty, col=hlines_col,
                   mgp=mgp.y)
            axis(side=2, at=hlines, mgp=mgp.y, tick=FALSE, las=las)

            if(!omit_points)
                points(match(geno, ugeno) + jitter, pheno, pch=pch, bg=bg, ...)

            if(seg_width > 0) {
                segments(x-seg_width/2, me, x+seg_width/2, me,
                         lwd=seg_lwd, col=seg_col)

                if(!is.null(SEmult)) {
                    segments(x, me[ugeno]-se[ugeno]*SEmult, x, me[ugeno]+se[ugeno]*SEmult,
                             lwd=seg_lwd, col=seg_col)
                    segments(x-seg_width/4, me-se*SEmult, x+seg_width/4, me-se*SEmult,
                             lwd=seg_lwd, col=seg_col)
                    segments(x-seg_width/4, me+se*SEmult, x+seg_width/4, me+se*SEmult,
                             lwd=seg_lwd, col=seg_col)
                }

            }

            box()
        }

    plot_pxg_internal(geno, pheno, bgcolor=bgcolor,
                      seg_width=seg_width, seg_lwd=seg_lwd, seg_col=seg_col,
                      hlines=hlines, hlines_col=hlines_col, hlines_lty=hlines_lty, hlines_lwd=hlines_lwd,
                      vlines_col=vlines_col, vlines_lty=vlines_lty, vlines_lwd=vlines_lwd, ...)

    # return mean and SE
    invisible(cbind(mean=me, SE=se))

}
