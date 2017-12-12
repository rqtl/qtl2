#' Get x-axis position for genomic location
#'
#' For a plot of [scan1()] results, get the x-axis
#' location that corresponds to a particular genomic location
#' (chromosome ID and position).
#'
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#' @param chr Selected chromosomes that were plotted (if used in the
#' call to [plot_scan1()].
#' @param gap The gap between chromosomes used in the call to
#' [plot_scan1()].
#' @param thechr Vector of chromosome IDs
#' @param thepos Vector of chromosomal positions
#'
#' @return A vector of x-axis locations.
#'
#' @details `thechr` and `thepos` should be the same length,
#' or should have length 1 (in which case they are expanded to the
#' length of the other vector).
#'
#' @export
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
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # plot the results for selected chromosomes
#' ylim <- c(0, maxlod(out)*1.02) # need to strip class to get overall max LOD
#' chr <- c(2,7,8,9,15,16)
#' plot(out, map, chr=chr, ylim=ylim)
#' plot(out, map, lodcolumn=2, chr=chr, col="violetred", add=TRUE)
#' legend("topleft", lwd=2, col=c("darkslateblue", "violetred"), colnames(out),
#'        bg="gray90")
#'
#' # Use xpos_scan1 to add points at the peaks
#' # first find the peaks with LOD > 3
#' peaks <- find_peaks(out, map)
#'
#' # keep just the peaks for chromosomes that were plotted
#' peaks <- peaks[peaks$chr %in% chr,]
#'
#' # find x-axis positions
#' xpos <- xpos_scan1(map, chr=chr, thechr=peaks$chr, thepos=peaks$pos)
#'
#' # point colors
#' ptcolor <- c("darkslateblue", "violetred")[match(peaks$lodcolumn, c("liver", "spleen"))]
#'
#' # plot points
#' points(xpos, peaks$lod, pch=21, bg=ptcolor)
xpos_scan1 <-
function(map, chr=NULL, gap=25, thechr, thepos)
{
    if(is.null(map)) stop("map is NULL")
    if(!is_nonneg_number(gap)) stop("gap should be a single non-negative number")

    # subset chromosomes
    if(!is.null(chr)) map <- map[chr]

    # start position for each chromosome
    start <- map_to_boundaries(map, gap)[1,] - vapply(map, "[", 1, 1)
    # if more than one chromosome, add gap/2
    if(length(map) > 1) start <- start + gap/2
    # add chromosome names
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
