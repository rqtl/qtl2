#' Plot comparison of two sets of genotype probabilities
#'
#' Plot a comparison of two sets of genotype probabilities for one individual on one chromosome, as a heat map.
#'
#' @param probs1 Genotype probabilities (as produced by \code{\link[qtl2geno]{calc_genoprob}})
#' or allele dosages (as produced by \code{\link[qtl2geno]{genoprob_to_alleleprob}}).
#' @param probs2 A second set of genotype probabilities, just like \code{probs1}.
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosome to plot; a single character string.
#' @param geno Optional vector of genotypes or alleles to be shown
#' (vector of integers or character strings)
#' @param threshold Threshold for genotype probabilities; only genotypes that achieve
#' this value somewhere on the chromosome (in one or the other set of probabilities) will be shown.
#' @param n_colors Number of colors in each color scale.
#' @param swap_axes If TRUE, swap the axes, so that the genotypes are
#' on the x-axis and the chromosome position is on the y-axis.
#' @param hlines Position of horizontal grid lines (use \code{NA} to avoid lines).
#' @param hlines_col Color of horizontal grid lines.
#' @param hlines_lwd Line width of horizontal grid lines.
#' @param hlines_lty Line type of horizontal grid lines.
#' @param vlines Position of vertical grid lines (use \code{NA} to avoid lines).
#' @param vlines_col Color of vertical grid lines.
#' @param vlines_lwd Line width of vertical grid lines.
#' @param vlines_lty Line type of vertical grid lines.
#' @param ... Additional graphics parameters passed to \code{\link[graphics]{image}}.
#'
#' @details
#' We plot the first set of probabilities in the range white to blue
#' and the second set in the range white to red and attempt to combine
#' them, for colors that are white, some amount of blue or red, or
#' where both are large something like blackish purple.
#'
#' @export
#' @importFrom graphics image par axis title box
#' @importFrom grDevices gray
plot_genoprobcomp <-
    function(probs1, probs2, map, ind=1, chr=NULL, geno=NULL,
             threshold=0, n_colors=256, swap_axes=FALSE,
             hlines=NULL, hlines_col="#B3B3B370", hlines_lwd=1, hlines_lty=1,
             vlines=NULL, vlines_col="#B3B3B370", vlines_lwd=1, vlines_lty=1,
             ...)
{
    # check inputs
    if(is.null(map)) stop("map is NULL")
    if(is.null(chr)) chr <- names(probs1)[1]
    if(length(chr) > 1) {
        warning("chr should have length 1; using the first value")
        chr <- chr[1]
    }
    if(!(chr %in% names(probs1) && chr %in% names(probs2)))
        stop("chr ", chr, " not found in both probs")
    if(!(chr %in% names(map))) stop("chr ", chr, " not found in map")
    if(length(ind) > 1) {
        warning("ind should have length 1; using the first value")
        ind <- ind[1]
    }

    # pull out selected chromosomes
    probs1 <- probs1[[chr]]
    probs2 <- probs2[[chr]]
    map <- map[[chr]]

    # find common markers
    mar1 <- dimnames(probs1)[[3]]
    mar2 <- dimnames(probs2)[[3]]
    marm <- names(map)
    mar <- mar1[mar1 %in% mar2 & mar1 %in% marm]
    if(length(mar) < 2) stop("<2 markers in common")
    probs1 <- probs1[,,mar,drop=FALSE]
    probs2 <- probs2[,,mar,drop=FALSE]
    map <- map[mar]

    # make sure they have the same genotypes
    if(ncol(probs1) != ncol(probs2) || !all(colnames(probs1) == colnames(probs2)))
        stop("Need ncol(probs1) == ncol(probs2) and same colnames")

    # pull out individual's probs; make it positions x probs
    if(is.character(ind) && !(ind %in% rownames(probs1) && ind %in% rownames(probs2)))
        stop("ind ", ind, " not found in both sets of probs")
    if(is.numeric(ind)) {
        if(nrow(probs1) != nrow(probs2))
            stop("If ind is numeric, must have nrow(probs1) == nrow(probs2)")
        if(ind < 1 || ind > nrow(probs1))
            stop("ind ", ind, " should be in the range [1, ", nrow(probs1), "]")
    }

    # now pull out individual and transpose
    probs1 <- t(probs1[ind,,])
    probs2 <- t(probs2[ind,,])

    # pull out selected genotypes
    if(!is.null(geno)) {
        if(is.numeric(geno)) {
            if(all(geno < 0)) {
                if(any(geno > -1 | geno < -ncol(probs1)))
                    stop("negative geno should be in the range [", -ncol(probs1), " , -1]")
                geno <- (1:ncol(probs1))[geno]
            }
            if(any(geno < 1 | geno > ncol(probs1)))
                stop("numeric geno should be in the range [1, ", ncol(probs1))
            geno <- colnames(probs1)[geno]
        }
        if(!all(geno %in% colnames(probs1)))
            stop("Not all geno in probs")
        probs1 <- probs1[,geno,drop=FALSE]
        probs2 <- probs2[,geno,drop=FALSE]
    }

    # drop genotypes that do exceed threshold
    if(threshold > 0) {
        geno_keep <- colSums(probs1 >= threshold) + colSums(probs2 >= threshold) > 0
        if(sum(geno_keep) == 0) stop("No genotype probabilities exceed the threshold")
        probs1 <- probs1[,geno_keep,drop=FALSE]
        probs2 <- probs2[,geno_keep,drop=FALSE]
    }

    # save dimnames
    dn <- dimnames(probs1)

    # split probabilities into integer levels
    #    need to screw around a bit with case of 1+epsilon
    probs1[probs1 > 1] <- 1
    probs1 <- matrix(as.numeric(cut(probs1, breaks=seq(0, 1, length=n_colors+1), include.lowest=TRUE, right=TRUE)),
                     ncol=ncol(probs1), nrow=nrow(probs1))

    probs2[probs2 > 1] <- 1
    probs2 <- matrix(as.numeric(cut(probs2, breaks=seq(0, 1, length=n_colors+1), include.lowest=TRUE, right=TRUE)),
                     ncol=ncol(probs2), nrow=nrow(probs2))

    # combine into one object
    probs <- (probs1-1)*n_colors + probs2
    dimnames(probs) <- dn

    # create color palette
    redmat <- bluemat <- matrix(ncol=3, nrow=n_colors)
    bluemat[,1] <- bluemat[,2] <- redmat[,2] <- redmat[,3] <- seq(1, 0, length=n_colors)
    bluemat[,3] <- redmat[,1] <- 1
    joint_colors <- 1:(n_colors^2)
    k <- 1
    for(i in 1:n_colors) {
        for(j in 1:n_colors) {
            joint_colors[k] <- rgb(bluemat[i,1]*redmat[j,1],
                                   bluemat[i,2]*redmat[j,2],
                                   bluemat[i,3]*redmat[j,3],
                                   maxColorValue=1)
            k <- k+1
        }
    }

    # separate positions if necessary
    tol <- 1e-6
    if(any(diff(map) < tol))
        map <- map + seq(0, tol, length.out=length(map))

    plot_genoprob_internal(probs, map, col=joint_colors, swap_axes=swap_axes,
                           hlines=hlines, hlines_col=hlines_col, hlines_lty=hlines_lty, hlines_lwd=hlines_lwd,
                           vlines=vlines, vlines_col=vlines_col, vlines_lty=vlines_lty, vlines_lwd=vlines_lwd,
                           zlim=c(1,n_colors^2),
                           ...)

}
