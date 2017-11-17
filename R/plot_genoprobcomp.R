#' Plot comparison of two sets of genotype probabilities
#'
#' Plot a comparison of two sets of genotype probabilities for one individual on one chromosome, as a heat map.
#'
#' @md
#'
#' @param probs1 Genotype probabilities (as produced by [qtl2geno::calc_genoprob()])
#' or allele dosages (as produced by [qtl2geno::genoprob_to_alleleprob()]).
#' @param probs2 A second set of genotype probabilities, just like `probs1`.
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
#' @param ... Additional graphics parameters passed to [graphics::image()].
#'
#' @details
#' We plot the first set of probabilities in the range white to blue
#' and the second set in the range white to red and attempt to combine
#' them, for colors that are white, some amount of blue or red, or
#' where both are large something like blackish purple.
#'
#' @seealso [plot_genoprob()]
#'
#' @examples
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' iron <- iron[1,"2"]   # subset to first individual on chr 2
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # in presence of a genotyping error, how much does error_prob matter?
#' iron$geno[[1]][1,3] <- 3
#' pr_e <- calc_genoprob(iron, map, error_prob=0.002)
#' pr_ne <- calc_genoprob(iron, map, error_prob=1e-15)
#'
#' # image of probabilities + comparison
#' par(mfrow=c(3,1))
#' plot_genoprob(pr_e, map, main="Allow errors")
#' plot_genoprob(pr_ne, map, main="Assume very low error error")
#' plot_genoprobcomp(pr_e, pr_ne, map, main="Comparison")
#'
#' @export
#' @importFrom grDevices rgb
plot_genoprobcomp <-
    function(probs1, probs2, map, ind=1, chr=NULL, geno=NULL,
             threshold=0, n_colors=256, swap_axes=FALSE, ...)
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
                geno <- seq_len(ncol(probs1))[geno]
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
    red <- blue <- matrix(1, ncol=3, nrow=n_colors^2)
    # using 0.2 - 1 for the range of colors so that the combination is more purple rather than black
    # (also if ending at 0, the blue is a bit ugly)
    blue[,1] <- blue[,2] <- rep(seq(1, 0.2, length=n_colors), each=n_colors)
    red[,2] <- red[,3] <-   rep(seq(1, 0.2, length=n_colors), n_colors)
    joint_colors <- apply(blue*red, 1, function(a) rgb(a[1], a[2], a[3], maxColorValue=1))

    plot_genoprob_internal(probs, map, col=joint_colors, swap_axes=swap_axes,
                           zlim=c(1,n_colors^2), ...)

}
