#' Plot phenotype vs genotype
#'
#' Plot phenotype vs genotype for a single putative QTL and a single phenotype.
#'
#' @param geno Vector of genotypes, for example as produced by
#' [maxmarg()] with specific `chr` and `pos`.
#' @param pheno Vector of phenotypes.
#' @param sort If TRUE, sort genotypes from largest to smallest.
#' @param SEmult If specified, interval estimates of the within-group
#' averages will be displayed, as `mean +/- SE * SEmult`.
#' @param pooledSD If TRUE and `SEmult` is specified, calculated
#' a pooled within-group SD. Otherwise, get separate estimates of
#' the within-group SD for each group.
#' @param swap_axes If TRUE, swap the axes, so that the genotypes are
#' on the y-axis and the phenotype is on the x-axis.
#' @param jitter Amount to jitter the points horizontally, if a vector
#' of length > 0, it is taken to be the actual jitter amounts
#' (with values between -0.5 and 0.5).
#' @param force_labels If TRUE, force all genotype labels to be shown.
#' @param alternate_labels If TRUE, place genotype labels in two rows
#' @param omit_points If TRUE, omit the points, just plotting the averages (and, potentially, the +/- SE intervals).
#' @param ... Additional graphics parameters, passed to `plot()`.
#'
#' @return (Invisibly) A matrix with rows being the genotype groups
#' and columns for the means and (if `SEmult` is specified) the SEs.
#'
#' @export
#' @importFrom graphics par segments title axis points
#' @importFrom stats lm runif sd
#'
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color, and
#' `seg_width`, `seg_lwd`, and `seg_col` to control the lines at the
#' confidence intervals. Further, `hlines`, `hlines_col`,
#' `hlines_lwd`, and `hlines_lty` to control the horizontal grid
#' lines. (Use `hlines=NA` to avoid plotting horizontal grid lines.)
#' Similarly `vlines`, `vlines_col`, `vlines_lwd`, and `vlines_lty`
#' for vertical grid lines. These are not included as formal
#' parameters in order to avoid cluttering the function definition.
#'
#' @seealso [plot_coef()]
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,"16"]}
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
#'
#' # include +/- 2 SE intervals
#' plot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          SEmult=2)
#'
#' # plot just the means
#' plot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          omit_points=TRUE)
#'
#' # plot just the means +/- 2 SEs
#' plot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          omit_points=TRUE, SEmult=2)
plot_pxg <-
    function(geno, pheno, sort=TRUE, SEmult=NULL, pooledSD=TRUE,
             swap_axes=FALSE, jitter=0.2,
             force_labels=TRUE, alternate_labels=FALSE,
             omit_points=FALSE, ...)
{
    if(is.null(geno)) stop("geno is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

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

    # force pheno to be numeric
    if(is.matrix(pheno) || is.data.frame(pheno)) {
        if(ncol(pheno)>1) {
            pheno <- pheno[,1,drop=FALSE]
            warning("pheno is a ", class(pheno)[1], "; using the first column")
        }
        pheno <- setNames(as.numeric(unlist(pheno)), rownames(pheno))
    }

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
        function(geno, pheno, swap_axes=FALSE, bgcolor="gray90",
                 seg_width=NULL, seg_lwd=2, seg_col="black",
                 hlines=NULL, hlines_col=NULL, hlines_lty=NULL, hlines_lwd=NULL,
                 vlines=NULL, vlines_col=NULL, vlines_lty=NULL, vlines_lwd=NULL,
                 xlim=NULL, ylim=NULL,
                 xaxs=NULL, yaxs=NULL, xlab=NULL, ylab=NULL,
                 mgp=c(2.6, 0.3, 0), mgp.x=mgp, mgp.y=mgp, las=NULL,
                 pch=21, bg="lightblue", main="", sub="", ...)
        {
            if(swap_axes) {
                if(is.null(ylim)) ylim <- c(0.5, length(ugeno)+0.5)

                # get x-axis limits
                if(is.null(xlim)) {
                    if(omit_points) {
                        xlim <- range(me, na.rm=TRUE)
                    } else {
                        xlim <- range(pheno, na.rm=TRUE)
                    }
                    if(!is.null(SEmult)) {
                        xlim <- range(c(xlim, me-se*SEmult, me+se*SEmult), na.rm=TRUE)
                    }
                }

                if(is.null(xlab)) xlab <- "Phenotype"
                if(is.null(ylab)) ylab <- "Genotype"

                if(is.null(yaxs)) yaxs <- "i"
                if(is.null(xaxs)) xaxs <- "r"

                if(is.null(las)) las <- 1

                if(is.null(hlines)) hlines <- seq_len(length(ugeno))
                if(is.null(vlines)) vlines <- pretty(pheno)
                if(is.null(hlines_col)) hlines_col <- "gray80"
                if(is.null(hlines_lwd)) hlines_lwd <- 3
                if(is.null(hlines_lty)) hlines_lty <- 1
                if(is.null(vlines_col)) vlines_col <- "white"
                if(is.null(vlines_lwd)) vlines_lwd <- 1
                if(is.null(vlines_lty)) vlines_lty <- 1

                x <- pheno
                y <- match(geno, ugeno) + jitter

            } else {
                if(is.null(xlim)) xlim <- c(0.5, length(ugeno)+0.5)

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

                if(is.null(ylab)) ylab <- "Phenotype"
                if(is.null(xlab)) xlab <- "Genotype"

                if(is.null(xaxs)) xaxs <- "i"
                if(is.null(yaxs)) yaxs <- "r"

                if(is.null(las)) {
                    las <- ifelse(length(ugeno) > 4, 2, 1)
                }

                if(is.null(hlines)) hlines <- pretty(pheno)
                if(is.null(vlines)) vlines <- seq_len(length(ugeno))
                if(is.null(vlines_col)) vlines_col <- "gray80"
                if(is.null(vlines_lwd)) vlines_lwd <- 3
                if(is.null(vlines_lty)) vlines_lty <- 1
                if(is.null(hlines_col)) hlines_col <- "white"
                if(is.null(hlines_lwd)) hlines_lwd <- 1
                if(is.null(hlines_lty)) hlines_lty <- 1

                y <- pheno
                x <- match(geno, ugeno) + jitter

            }

            plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n",
                 yaxt="n", xaxs=xaxs, yaxs=yaxs, main=main, sub=sub)
            title(xlab=xlab, mgp=mgp.x)
            title(ylab=ylab, mgp=mgp.y)

            u <- par("usr")
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)

            # grid lines (if swap_axes, horizontal over vertical
            if(!swap_axes) {
                if(!(length(vlines)==1 && is.na(vlines)))
                    abline(v=vlines, lwd=vlines_lwd, lty=vlines_lty, col=vlines_col)
            }
            if(!(length(hlines)==1 && is.na(hlines)))
                abline(h=hlines, lwd=hlines_lwd, lty=hlines_lty, col=hlines_col)
            if(swap_axes) {
                if(!(length(vlines)==1 && is.na(vlines)))
                    abline(v=vlines, lwd=vlines_lwd, lty=vlines_lty, col=vlines_col)
            }

            # x-axis labels
            xtick <- seq_along(ugeno)
            side <- swap_axes + 1
            this_mgp <- list(mgp.x, mgp.y)[[side]]
            if(force_labels) {
                    if(alternate_labels) {
                        for(i in xtick)
                            axis(side=side, at=i, labels=ugeno[i], mgp=this_mgp, tick=FALSE,
                                 line=2 - (i %% 2)-1, las=las)
                    } else {
                        for(i in xtick)
                            axis(side=side, at=i, labels=ugeno[i], mgp=this_mgp, tick=FALSE, las=las)
                    }
            } else {
                if(alternate_labels) {
                    odd <- seq(1, length(ugeno), by=2)
                    axis(side=side, at=xtick[odd], labels=ugeno[odd], mgp=this_mgp, tick=FALSE, las=las)
                    if(length(ugeno) > 1) {
                        axis(side=side, at=xtick[-odd], labels=ugeno[-odd], mgp=this_mgp, tick=FALSE,
                             line=1, las=las)
                    }
                } else {
                    axis(side=side, at=xtick, labels=ugeno, mgp=this_mgp, tick=FALSE, las=las)
                }
            }

            # add points
            if(!omit_points)
                points(x, y, pch=pch, bg=bg, ...)

            # y-axis labels
            if(swap_axes) {
                axis(side=1, at=vlines, mgp=mgp.x, tick=FALSE, las=las)
            } else {
                axis(side=2, at=hlines, mgp=mgp.y, tick=FALSE, las=las)
            }

            # determine default CI segment width; use width of plot but make it no bigger than 0.2 and no smaller than 0.05
            if(is.null(seg_width)) {
                seg_width <- ifelse(swap_axes, diff(par("usr"))[3:4], diff(par("usr"))[1:2])/20
                if(seg_width > 0.2) seg_width <- 0.2
                if(seg_width < 0.05) seg_width <- 0.05
            }

            # add CIs
            if(seg_width > 0) {
                if(swap_axes) {
                    segments(me, xtick-seg_width/2, me, xtick+seg_width/2,
                             lwd=seg_lwd, col=seg_col)
                }
                else {
                    segments(xtick-seg_width/2, me, xtick+seg_width/2, me,
                             lwd=seg_lwd, col=seg_col)
                }

                if(!is.null(SEmult)) {
                    if(swap_axes) {
                        segments(me[ugeno]+se[ugeno]*SEmult, xtick, me[ugeno]-se[ugeno]*SEmult, xtick,
                                 lwd=seg_lwd, col=seg_col)
                    } else {
                        segments(xtick, me[ugeno]-se[ugeno]*SEmult, xtick, me[ugeno]+se[ugeno]*SEmult,
                                 lwd=seg_lwd, col=seg_col)
                    }
                }

            }

            box()
        }

    plot_pxg_internal(geno, pheno, swap_axes=swap_axes, ...)


    # return mean and SE
    invisible(cbind(mean=me, SE=se))

}
