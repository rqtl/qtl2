# find_peaks_and_lodint
#
# This is like find_peaks but the results also include lod support intervals.
# The same input as find_peaks, which calls this function when drop is given.
find_peaks_and_lodint <-
    function(scan1_output, threshold=3, peakdrop=Inf, drop=1.8,
             thresholdX=NULL, peakdropX=NULL, dropX=NULL, cores=1)
{
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

    # make the drops have length n_lod
    if(length(drop)==1) drop <- rep(drop, n_lod)
    if(length(drop) != n_lod)
        stop("drop should have length 1 or ", n_lod)
    if(is.null(dropX)) dropX <- drop
    if(length(dropX)==1) dropX <- rep(dropX, n_lod)
    if(length(dropX) != n_lod)
        stop("dropX should have length 1 or ", n_lod)

    if(any(drop > peakdrop | dropX > peakdropX))
        stop("Must have drop <= peakdrop")

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

        peaks <- .find_peaks_and_lodint(this_lod,
                                        ifelse(is_x_chr[chr], thresholdX[lodcol], threshold[lodcol]),
                                        ifelse(is_x_chr[chr], peakdropX[lodcol], peakdrop[lodcol]),
                                        ifelse(is_x_chr[chr], dropX[lodcol], drop[lodcol]))

        n_peaks <- length(peaks)
        if(n_peaks==0) return(NULL)

        # calculate peak positions
        # (this deals with ties in the LOD scores:
        #  take average of multiple positions sharing the maximum LOD score)
        # and remember that .find_peaks returns indexes starting at 0
        ci_lo <- vapply(peaks, function(a) map[[chr]][a[1]+1], 0.0)
        ci_hi <- vapply(peaks, function(a) map[[chr]][a[2]+1], 0.0)
        peak_pos <- vapply(peaks, function(a) mean(map[[chr]][a[-(1:2)]+1]), 0.0)
        peak_lod <- vapply(peaks, function(a) this_lod[a[3]+1], 0.0)

        data.frame(lodindex=rep(lodcol,n_peaks),
                   lodcolumn=rep(lodnames[lodcol],n_peaks),
                   chr=rep(names(map)[chr],n_peaks),
                   pos=peak_pos,
                   lod=peak_lod,
                   ci_lo=ci_lo,
                   ci_hi=ci_hi,
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
