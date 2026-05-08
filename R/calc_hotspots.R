#' Calculate QTL hotspots
#'
#' For a set of QTL locations, calculate a running count in a sliding
#' window across the genome.
#'
#' @param peaks Data frame of QTL results, as output by [find_peaks()]
#'     Should contain columns `chr` and `pos`.
#'
#' @param map Marker map, as a list of chromosomes, each being a vector of
#'     positions. Hotspot counts will be calculated at these
#'     positions.
#'
#' @param window Window size for counting QTL.
#'
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @param quiet If TRUE, don't print any messages.
#'
#' @return An object of class `"scan1"`: a matrix with a single
#'     column, of counts, with rownames being the marker names in
#'     `map`. The column name is `"num_qtl"`.
#'
#'
#' @seealso [find_peaks()], [plot_lodpeaks()], [plot_cistrans()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example pQTL results (from Keele et al. 2026, https://doi.org/10.1016/j.xgen.2025.101069)
#' # contains qtl, map_endpoints, and pheno_pos
#' url <- "https://kbroman.org/qtl2/assets/sampledata/pqtl_data.RData"
#' tempfile <- file.path(tempdir(), basename(url))
#' download.file(url, tempfile)
#' load(tempfile)
#' unlink(tempfile)
#'
#' hotspots <- calc_hotspots(qtl, map, window=2)
#' plot(hotspots, map, ylab="No. QTL")
#' find_peaks(hotspots, map, threshold=20)
#' }

calc_hotspots <-
    function(peaks, map, window=1, cores=1, quiet=TRUE)
{
    if(!is.data.frame(peaks) || !all(c("chr", "pos") %in% names(peaks))) {
        stop("peaks should be a data frame with columns chr and pos")
    }
    stopifnot(is.list(map))

    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    bychr_func <-
        function(chr)
    {
        result <- .running_count(peaks$pos[peaks$chr==chr], map[[chr]], window)

        names(result) <- names(map[[chr]])
        result
    }

    result <- as.matrix(unlist(cluster_lapply(cores, names(map), bychr_func)))

    class(result) <- c("scan1", "matrix")
    colnames(result) <- "num_qtl"
    result
}
