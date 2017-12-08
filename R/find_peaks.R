# find_peaks
#' Find peaks in a set of LOD curves
#'
#' Find peaks in a set of LOD curves (output from [scan1()]
#'
#' @md
#'
#' @param scan1_output An object of class `"scan1"` as returned by
#' [scan1()].
#' @param map A list of vectors of marker positions, as produced by
#' [insert_pseudomarkers()].
#' @param threshold Minimum LOD score for a peak (can be a vector with
#' separate thresholds for each lod score column in
#' `scan1_output`)
#' @param peakdrop Amount that the LOD score must drop between peaks,
#' if multiple peaks are to be defined on a chromosome. (Can be a vector with
#' separate values for each lod score column in
#' `scan1_output`.)
#' @param drop If provided, LOD support intervals are included in the
#' results, and this indicates the amount to drop in the support
#' interval. (Can be a vector with
#' separate values for each lod score column in
#' `scan1_output`.) Must be \eqn{\le} `peakdrop`
#' @param prob If provided, Bayes credible intervals are included in the
#' results, and this indicates the nominal coverage.
#' (Can be a vector with
#' separate values for each lod score column in
#' `scan1_output`.) Provide just one of `drop` and `prob`.
#' @param thresholdX Separate threshold for the X chromosome; if
#' unspecified, the same threshold is used for both autosomes and the
#' X chromosome. (Like `threshold`, this can be a vector with
#' separate thresholds for each lod score column.)
#' @param peakdropX Like `peakdrop`, but for the X chromosome; if
#' unspecified, the same value is used for both autosomes and the X
#' chromosome.  (Can be a vector with separate values for each lod
#' score column in `scan1_output`.)
#' @param dropX Amount to drop for LOD support intervals on the X
#' chromosome.  Ignored if `drop` is not provided. (Can be a
#' vector with separate values for each lod score column in
#' `scan1_output`.)
#' @param probX Nominal coverage for Bayes intervals on the X
#' chromosome.  Ignored if `prob` is not provided. (Can be a
#' vector with separate values for each lod score column in
#' `scan1_output`.)
#' @param expand2markers If TRUE (and if `drop` or `prob` is
#' provided, so that QTL intervals are calculated), QTL intervals are
#' expanded so that their endpoints are at genetic markers.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A data frame with each row being a single peak on a single
#' chromosome for a single LOD score column, and with columns
#' * `lodindex` - lod column index
#' * `lodcolumn` - lod column name
#' * `chr` - chromosome ID
#' * `pos` - peak position
#' * `lod` - lod score at peak
#'
#' If `drop` or `prob` is provided, the results will include
#' two additional columns: `ci_lo` and `ci_hi`, with the
#' endpoints of the LOD support intervals or Bayes credible wintervals.
#'
#' @details For each lod score column on each chromosome, we return a
#' set of peaks defined as local maxima that exceed the specified
#' `threshold`, with the requirement that the LOD score must have
#' dropped by at least `peakdrop` below the lowest of any two
#' adjacent peaks.
#'
#' At a given peak, if there are ties, with multiple positions jointly
#' achieving the maximum LOD score, we take the average of these
#' positions as the location of the peak.
#'
#' @export
#'
#' @seealso [scan1()], [lod_int()], [bayes_int()]
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
#' # find just the highest peak on each chromosome
#' find_peaks(out, map, threshold=3, peakdrop=Inf)
#'
#' # possibly multiple peaks per chromosome
#' find_peaks(out, map, threshold=3, peakdrop=2)
#'
#' # possibly multiple peaks, also getting LOD support intervals
#' find_peaks(out, map, threshold=3, peakdrop=2, drop=1.5)
#'
#' # possibly multiple peaks, also getting Bayes intervals
#' find_peaks(out, map, threshold=3, peakdrop=2, prob=0.95)
find_peaks <-
    function(scan1_output, map, threshold=3, peakdrop=Inf, drop=NULL, prob=NULL,
             thresholdX=NULL, peakdropX=NULL, dropX=NULL, probX=NULL,
             expand2markers=TRUE, cores=1)
{
    # align scan1_output and map
    tmp <- align_scan1_map(scan1_output, map)
    scan1_output <- tmp$scan1
    map <- tmp$map

    if(nrow(scan1_output) != length(unlist(map)))
        stop("nrow(scan1_output) [", nrow(scan1_output), "] != number of positions in map [",
             length(unlist(map)), "]")

    if(!is.null(drop) && !is.null(prob))
        stop('No more than one of "drop" and "prob" should be provided')

    if(!is.null(drop)) # also include lod support intervals
        return(find_peaks_and_lodint(scan1_output, map, threshold, peakdrop, drop,
                                     thresholdX, peakdropX, dropX,
                                     expand2markers, cores))
    if(!is.null(prob)) # also include Bayes credible intervals
        return(find_peaks_and_bayesint(scan1_output, map, threshold, peakdrop, prob,
                                       thresholdX, peakdropX, probX,
                                       expand2markers, cores))

    if(!is.null(dropX))
        warning("dropX ignored if drop is not provided")
    if(!is.null(probX))
        warning("probX ignored if prob is not provided")

    lodnames <- colnames(scan1_output)
    n_lod <- length(lodnames)

    # make the thresholds have length n_lod
    if(length(threshold)==1) threshold <- rep(threshold, n_lod)
    if(length(threshold) != n_lod)
        stop("threshold should have length 1 or ", n_lod)
    if(is.null(thresholdX)) thresholdX <- threshold
    if(length(thresholdX)==1) thresholdX <- rep(thresholdX, n_lod)
    if(length(thresholdX) != n_lod)
        stop("thresholdX should have length 1 or ", n_lod)

    # make the peakdrops have length n_lod
    if(length(peakdrop)==1) peakdrop <- rep(peakdrop, n_lod)
    if(length(peakdrop) != n_lod)
        stop("peakdrop should have length 1 or ", n_lod)
    if(is.null(peakdropX)) peakdropX <- peakdrop
    if(length(peakdropX)==1) peakdropX <- rep(peakdropX, n_lod)
    if(length(peakdropX) != n_lod)
        stop("peakdropX should have length 1 or ", n_lod)

    # split lod scores by column and by chromosome
    lod <- as.list(as.data.frame(unclass(scan1_output)))
    chr <- rep(seq_along(map), vapply(map, length, 1))
    lod <- lapply(lod, split, chr)

    # batch info for parallel-processing
    batch <- cbind(rep(seq_along(lod), vapply(lod, length, 1)),
                   unlist(lapply(lod, function(a) seq_along(a))))
    dimnames(batch) <- NULL
    batch_index <- seq_len(nrow(batch))

    # set up parallel analysis
    cores <- setup_cluster(cores)

    # X chr info
    is_x_chr <- attr(map, "is_x_chr")
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(map))

    # function to be applied to each batch
    by_batch_func <- function(i) {
        lodcol <- batch[i,1]
        chr <- batch[i,2]

        this_lod <- lod[[lodcol]][[chr]]

        peaks <- .find_peaks(this_lod,
                             ifelse(is_x_chr[chr], thresholdX[lodcol], threshold[lodcol]),
                             ifelse(is_x_chr[chr], peakdropX[lodcol], peakdrop[lodcol]))

        n_peaks <- length(peaks)
        if(n_peaks==0) return(NULL)

        # calculate peak positions
        # (this deals with ties in the LOD scores:
        #  take average of multiple positions sharing the maximum LOD score)
        # and remember that .find_peaks returns indexes starting at 0
        peak_pos <- vapply(peaks, function(a) mean(map[[chr]][a+1]), 0.0)
        peak_lod <- vapply(peaks, function(a) this_lod[a[1]+1], 0.0)

        data.frame(lodindex=rep(lodcol,n_peaks),
                   lodcolumn=rep(lodnames[lodcol],n_peaks),
                   chr=rep(names(map)[chr],n_peaks),
                   pos=peak_pos,
                   lod=peak_lod,
                   row.names=seq_along(peaks),
                   stringsAsFactors=FALSE)
    }

    if(n_cores(cores)==1) {
        peaks <- lapply(batch_index, by_batch_func)
    } else {
        peaks <- cluster_lapply(cores, batch_index, by_batch_func)
    }

    result <- NULL
    for(p in peaks) result <- rbind(result, p)

    rownames(result) <- NULL
    result$chr <- factor(result$chr, names(map))
    result
}
