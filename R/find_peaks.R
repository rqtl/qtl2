# find_peaks
#' Find peaks in a set of LOD curves
#'
#' Find peaks in a set of LOD curves (output from \code{\link{scan1}}
#'
#' @param scan1_output An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}}.
#' @param threshold Minimum LOD score for a peak (can be a vector with
#' separate thresholds for each lod score column in
#' \code{scan1_output})
#' @param peakdrop Amount that the LOD score must drop between peaks,
#' if multiple peaks are to be defined on a chromosome. (Can be a vector with
#' separate values for each lod score column in
#' \code{scan1_output}.)
#' @param drop If provided, LOD support intervals are included in the
#' results, and this indicates the amount to drop in the support
#' interval. (Can be a vector with
#' separate values for each lod score column in
#' \code{scan1_output}.) Must be \eqn{\le} \code{peakdrop}
#' @param thresholdX Separate threshold for the X chromosome; if
#' unspecified, the same threshold is used for both autosomes and the
#' X chromosome. (Like \code{threshold}, this can be a vector with
#' separate thresholds for each lod score column.)
#' @param peakdropX Like \code{peakdrop}, but for the X chromosome; if
#' unspecified, the same value is used for both autosomes and the X
#' chromosome.  (Can be a vector with separate values for each lod
#' score column in \code{scan1_output}.)
#' @param dropX Amount to drop for LOD support intervals on the X
#' chromosome.  Ignored if \code{drop} is not provided. (Can be a
#' vector with separate values for each lod score column in
#' \code{scan1_output}.)
#' @param expand2markers If TRUE (and if \code{drop} is provided, so
#' that QTL intervals are calculated), QTL intervals are expanded so
#' that their endpoints are at genetic markers.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A data frame with each row being a single peak on a single
#' chromosome for a single LOD score column, and with columns
#' \itemize{
#' \item \code{lodindex} - lod column index
#' \item \code{lodcolumn} - lod column name
#' \item \code{chr} - chromosome ID
#' \item \code{pos} - peak position
#' \item \code{lod} - lod score at peak
#' }
#'
#' If \code{drop} is provided, the results will include two additional
#' columns: \code{ci_lo} and \code{ci_hi}, with the endpoints of the
#' LOD support intervals.
#'
#' @details For each lod score column on each chromosome, we return a
#' set of peaks defined as local maxima that exceed the specified
#' \code{threshold}, with the requirement that the LOD score must have
#' dropped by at least \code{peakdrop} below the lowest of any two
#' adjacent peaks.
#'
#' At a given peak, if there are ties, with multiple positions jointly
#' achieving the maximum LOD score, we take the average of these
#' positions as the location of the peak.
#'
#' @export
#'
#' @seealso \code{\link{scan1}}
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
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # find just the highest peak on each chromosome
#' find_peaks(out, threshold=3, peakdrop=Inf)
#'
#' # possibly multiple peaks per chromosome
#' find_peaks(out, threshold=3, peakdrop=2)
#'
#' # possibly multiple peaks, also getting LOD support intervals
#' find_peaks(out, threshold=3, peakdrop=2, drop=1.5)
find_peaks <-
    function(scan1_output, threshold=3, peakdrop=Inf, drop=NULL,
             thresholdX=NULL, peakdropX=NULL, dropX=NULL,
             expand2markers=TRUE, cores=1)
{
    if(!is.null(drop)) # also include lod support intervals
        return(find_peaks_and_lodint(scan1_output, threshold, peakdrop, drop,
                                     thresholdX, peakdropX, dropX,
                                     expand2markers, cores))

    lodnames <- colnames(scan1_output$lod)
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
    lod <- as.list(as.data.frame(scan1_output$lod))
    map <- scan1_output$map
    chr <- rep(seq(along=map), vapply(map, length, 1))
    lod <- lapply(lod, split, chr)

    # batch info for parallel-processing
    batch <- cbind(rep(seq(along=lod), vapply(lod, length, 1)),
                   unlist(lapply(lod, function(a) seq(along=a))))
    dimnames(batch) <- NULL
    batch_index <- 1:nrow(batch)

    # set up parallel analysis
    cores <- setup_cluster(cores)

    # X chr info
    is_x_chr <- scan1_output$is_x_chr

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
                   row.names=seq(along=peaks),
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
    result
}
