# bayes_int
#' Calculate Bayes credible intervals
#'
#' Calculate Bayes credible intervals for a single LOD curve on a single
#' chromosome, with the ability to identify intervals for multiple LOD
#' peaks.
#'
#' @md
#'
#' @param scan1_output An object of class `"scan1"` as returned by
#' [scan1()].
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#' @param chr Chromosome ID to consider (must be a single value).
#' @param lodcolumn LOD score column to consider (must be a single value).
#' @param threshold Minimum LOD score for a peak.
#' @param peakdrop Amount that the LOD score must drop between peaks,
#' if multiple peaks are to be defined on a chromosome.
#' @param prob Nominal coverage for the interval.
#' @param expand2markers If TRUE, QTL intervals are expanded so
#' that their endpoints are at genetic markers.
#'
#' @return A matrix with three columns:
#' * `ci_lo` - lower bound of interval
#' * `pos` - peak position
#' * `ci_hi` - upper bound of interval
#'
#' Each row corresponds to a different peak.
#'
#' @details We identify a set of peaks defined as local maxima that
#' exceed the specified `threshold`, with the requirement that
#' the LOD score must have dropped by at least `peakdrop` below
#' the lowest of any two adjacent peaks.
#'
#' At a given peak, if there are ties, with multiple positions jointly
#' achieving the maximum LOD score, we take the average of these
#' positions as the location of the peak.
#'
#' The default is to use `threshold=0`, `peakdrop=Inf`, and
#' `prob=0.95`. We then return results a single peak, no matter the
#' maximum LOD score, and give a 95% Bayes credible interval.
#'
#' @export
#'
#' @seealso [lod_int()], [find_peaks()], [scan1()]
#'
#' @examples
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,7]}
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
#' # 95% Bayes credible interval for QTL on chr 7, first phenotype
#' bayes_int(out, map, chr=7, lodcolum=1)
bayes_int <-
    function(scan1_output, map, chr, lodcolumn=1, threshold=0,
             peakdrop=Inf, prob=0.95, expand2markers=TRUE)
{
    # align scan1_output and map
    tmp <- align_scan1_map(scan1_output, map)
    scan1_output <- tmp$scan1
    map <- tmp$map

    if(nrow(scan1_output) != length(unlist(map)))
        stop("nrow(scan1_output) [", nrow(scan1_output), "] != number of positions in map [",
             length(unlist(map)), "]")

    if(missing(chr) || is.null(chr)) { # just use the first chr
        chr <- names(map)[1]
    }

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
    if(length(prob) > 1) {
        warning("prob should have length 1; using the first value")
        prob <- prob[1]
    }

    if(lodcolumn < 1 || lodcolumn > ncol(scan1_output))
        stop("lodcolumn should be between 1 and ", ncol(scan1_output))

    # for chr to be character string
    chr <- as.character(chr)

    if(!(chr %in% names(map)))
        stop("Chromosome ", chr, " not found")

    scan1_output <- subset_scan1(scan1_output, map, chr, lodcolumn)

    lod <- unclass(scan1_output)
    map <- map[[chr]]

    peaks <- .find_peaks_and_bayesint(lod, map, threshold, peakdrop, prob)
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
    for(i in seq_len(n_peaks)) {
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
    rownames(result) <- seq_len(nrow(result))

    result
}
