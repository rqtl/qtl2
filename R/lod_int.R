# lod_int
#' Calculate LOD support intervals
#'
#' Calculate LOD support intervals for a single LOD curve on a single
#' chromosome, with the ability to identify intervals for multiple LOD
#' peaks.
#'
#' @param scan1_output An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}}.
#' @param chr Chromosome ID to consider (must be a single value).
#' @param lodcolumn LOD score column to consider (must be a single value).
#' @param threshold Minimum LOD score for a peak.
#' @param peakdrop Amount that the LOD score must drop between peaks,
#' if multiple peaks are to be defined on a chromosome.
#' @param drop Amount to drop in the support interval.  Must be
#' \eqn{\le} \code{peakdrop}
#' @param expand2markers If TRUE, QTL intervals are expanded so
#' that their endpoints are at genetic markers.
#'
#' @return A matrix with three columns:
#' \itemize{
#' \item \code{ci_lo} - lower bound of interval
#' \item \code{pos} - peak position
#' \item \code{ci_hi} - upper bound of interval
#' }
#' Each row corresponds to a different peak.
#'
#' @details We identify a set of peaks defined as local maxima that
#' exceed the specified \code{threshold}, with the requirement that
#' the LOD score must have dropped by at least \code{peakdrop} below
#' the lowest of any two adjacent peaks.
#'
#' At a given peak, if there are ties, with multiple positions jointly
#' achieving the maximum LOD score, we take the average of these
#' positions as the location of the peak.
#'
#' The default is to use \code{threshold=0}, \code{peakdrop=Inf}, and
#' \code{drop=1.5}. We then return results a single peak, no matter the
#' maximum LOD score, and give a 1.5-LOD support interval.
#'
#' @export
#'
#' @seealso \code{\link{bayes_int}}, \code{\link{find_peaks}}, \code{\link{scan1}}
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' \dontshow{iron <- iron[,7]}
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
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # 1.5-LOD support interval for QTL on chr 7, first phenotype
#' lod_int(out, chr=7, lodcolum=1)
lod_int <-
    function(scan1_output, chr, lodcolumn=1, threshold=0,
             peakdrop=Inf, drop=1.5, expand2markers=TRUE)
{
    if(length(chr) > 1) {
        warning("chr should have length 1; using the first value")
        chr <- chr[1]
    }
    if(length(lodcolumn) > 1) {
        warning("lodcolumn should have length 1; using the first value")
        lodcolumn <- lodcolumn[1]
    }
    if(length(threshold) > 1) {
        warning("threshold should have length 1; using the first value")
        threshold <- threshold[1]
    }
    if(length(peakdrop) > 1) {
        warning("peakdrop should have length 1; using the first value")
        peakdrop <- peakdrop[1]
    }
    if(length(drop) > 1) {
        warning("drop should have length 1; using the first value")
        drop <- drop[1]
    }

    if(drop > peakdrop)
        stop("Must have drop <= peakdrop")

    scan1_output <- scan1_output[chr, lodcolumn]

    lod <- scan1_output$lod
    map <- scan1_output$map[[1]]

    peaks <- .find_peaks_and_lodint(lod, threshold, peakdrop, drop)
    n_peaks <- length(peaks)
    if(n_peaks==0) {
        result <- matrix(nrow=0, ncol=3)
        colnames(result) <- c("ci_lo", "pos", "ci_hi")
        return(result)
    }

    # calculate peak positions
    # (this deals with ties in the LOD scores:
    #  take average of multiple positions sharing the maximum LOD score)
    # and remember that .find_peaks returns indexes starting at 0
    ci_lo <- ci_hi <- rep(0, n_peaks)
    for(i in 1:n_peaks) {
        if(expand2markers)
            ci <- expand_interval_to_markers(peaks[[i]][1:2]+1, map)
        else
            ci <- map[peaks[[i]][1:2]+1]
        ci_lo[i] <- ci[1]
        ci_hi[i] <- ci[2]
    }
    peak_pos <- vapply(peaks, function(a) mean(map[a[-(1:2)]+1]), 0.0)

    result <- cbind(ci_lo=ci_lo,
                    pos=peak_pos,
                    ci_hi=ci_hi)
    rownames(result) <- 1:nrow(result)

    result
}
