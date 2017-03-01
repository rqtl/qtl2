# find_peaks_and_bayesint
#
# This is like find_peaks but the results also include Bayes credible intervals.
# The same input as find_peaks, which calls this function when prob is given.
find_peaks_and_bayesint <-
    function(scan1_output, map, threshold=3, peakdrop=Inf, prob=0.95,
             thresholdX=NULL, peakdropX=NULL, probX=NULL,
             expand2markers=TRUE, cores=1)
{
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

    # make the probs have length n_lod
    if(length(prob)==1) prob <- rep(prob, n_lod)
    if(length(prob) != n_lod)
        stop("prob should have length 1 or ", n_lod)
    if(is.null(probX)) probX <- prob
    if(length(probX)==1) probX <- rep(probX, n_lod)
    if(length(probX) != n_lod)
        stop("probX should have length 1 or ", n_lod)

    # split lod scores by column and by chromosome
    lod <- as.list(as.data.frame(unclass(scan1_output)))
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
    is_x_chr <- attr(map, "is_x_chr")
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(map))

    # function to be applied to each batch
    by_batch_func <- function(i) {
        lodcol <- batch[i,1]
        chr <- batch[i,2]

        this_lod <- lod[[lodcol]][[chr]]

        peaks <- .find_peaks_and_bayesint(this_lod, map[[chr]],
                                          ifelse(is_x_chr[chr], thresholdX[lodcol], threshold[lodcol]),
                                          ifelse(is_x_chr[chr], peakdropX[lodcol], peakdrop[lodcol]),
                                          ifelse(is_x_chr[chr], probX[lodcol], prob[lodcol]))

        n_peaks <- length(peaks)
        if(n_peaks==0) return(NULL)

        # calculate peak positions
        # (this deals with ties in the LOD scores:
        #  take average of multiple positions sharing the maximum LOD score)
        # and remember that .find_peaks returns indexes starting at 0
        ci_lo <- ci_hi <- rep(0, n_peaks)
        for(i in 1:n_peaks) {
            if(expand2markers)
                ci <- expand_interval_to_markers(peaks[[i]][1:2]+1, map[[chr]])
            else
                ci <- map[[chr]][peaks[[i]][1:2]+1]
            ci_lo[i] <- ci[1]
            ci_hi[i] <- ci[2]
        }
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
