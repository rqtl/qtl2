# find_peaks
#' Find peaks in a set of LOD curves
#'
#' Find peaks in a set of LOD curves (output from \code{\link{scan1}}
#'
#' @param scan1_output An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}}.
#' @param threshold Minimum LOD score for a peak
#' @param peakdrop Amount that the LOD score must drop between peaks,
#' if multiple peaks are to be defined on a chromosome.
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
#' @details For each lod score column on each chromosome, we return a
#' set of peaks defined as local maxima that exceed the specified
#' \code{threshold}, with the requirement that the LOD score must have
#' dropped by at least \code{peakdrop} below the lowest of any two
#' adjacent peaks.
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
find_peaks <-
    function(scan1_output, threshold=3, peakdrop=Inf,
             cores=1)
{
    # split lod scores by column and by chromosome
    lodnames <- colnames(scan1_output$lod)
    lod <- as.list(as.data.frame(scan1_output$lod))
    map <- scan1_output$map
    chr <- rep(seq(along=map), vapply(map, length, 1))
    lod <- lapply(lod, split, chr)

    batch <- cbind(rep(seq(along=lod), vapply(lod, length, 1)),
                   unlist(lapply(lod, function(a) seq(along=a))))
    dimnames(batch) <- NULL
    batch_index <- 1:nrow(batch)

    # set up parallel analysis
    cores <- setup_cluster(cores)

    by_batch_func <- function(i) {
        lodcoli <- batch[i,1]
        chri <- batch[i,2]

        this_lod <- lod[[lodcoli]][[chri]]

        peaks <- .find_peaks(this_lod, threshold, peakdrop)

        n_peaks <- length(peaks)
        if(n_peaks==0) return(NULL)

        data.frame(lodindex=rep(lodcoli,n_peaks),
                   lodcolumn=rep(lodnames[lodcoli],n_peaks),
                   chr=rep(names(map)[chri],n_peaks),
                   pos=map[[chri]][peaks],
                   lod=this_lod[peaks],
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
