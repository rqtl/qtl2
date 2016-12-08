#' Get x-axis position for genomic location
#'
#' For a plot of \code{\link[qtl2scan]{scan1}} results, get the x-axis
#' location that corresponds to a particular genomic location
#' (chromosome ID and position).
#'
#' @param scan1_output Output of \code{\link[qtl2scan]{scan1}} that
#' was used in the call to \code{\link{plot_scan1}}. This can also be
#' a map object (a list of vectors of positions).
#' @param chr Selected chromosomes that were plotted (if used in the
#' call to \code{\link{plot_scan1}}).
#' @param gap The gap between chromosomes used in the call to
#' \code{\link{plot_scan1}}.
#' @param thechr Vector of chromosome IDs
#' @param thepos Vector of chromosomal positions
#'
#' @return A vector of x-axis locations.
#'
#' @details \code{thechr} and \code{thepos} should be the same length,
#' or should have length 1 (in which case they are expanded to the
#' length of the other vector).
#'
#' @export
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
#' # Use xpos_scan1 to add points at the peaks
#' # first find the peaks with LOD > 3
#' peaks <- find_peaks(out)
#'
#' # keep just the peaks for chromosomes that were plotted
#' peaks <- peaks[peaks$chr %in% chr,]
#'
#' # find x-axis positions
#' xpos <- xpos_scan1(out, chr=chr, thechr=peaks$chr, thepos=peaks$pos)
#'
#' # point colors
#' ptcolor <- c("darkslateblue", "violetred")[match(peaks$lodcolumn, c("liver", "spleen"))]
#'
#' # plot points
#' points(xpos, peaks$lod, pch=21, bg=ptcolor)
xpos_scan1 <-
function(scan1_output, chr=NULL, gap=25, thechr, thepos)
{
    if("map" %in% names(scan1_output))
        map <- scan1_output$map
    else map <- scan1_output

    # subset chromosomes
    if(!is.null(chr)) map <- map[chr]

    # start position for each chromosome (need to add gap/2 and subtract off the initial chromosome positions)
    start <- map_to_boundaries(map, gap)[1,] - vapply(map, "[", 1, 1) + gap/2
    names(start) <- names(map)

    if(length(thechr)==1)
        thechr <- rep(thechr, length(thepos))
    if(length(thepos)==1)
        thepos <- rep(thepos, length(thechr))
    if(length(thechr) != length(thepos))
        stop("thechr (length ", length(thechr), ") and thepos (length ",
             length(thepos), ") must be the same length, or length 1")

    uchr <- unique(thechr)
    if(!all(uchr %in% names(start)))
        stop("Unknown chromosome IDs: ",
             paste(uchr[!(uchr %in% names(start))], collapse=", "))

    result <- start[thechr] + thepos
    names(result) <- NULL

    result
}
