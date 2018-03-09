#' Locate crossovers
#'
#' Estimate the locations of crossovers in each individual on each chromosome.
#'
#' @md
#'
#' @param geno List of matrices of genotypes (output of [maxmarg()] or [viterbi()]).
#' @param map List of vectors with the map positions of the markers.
#' @param quiet If FALSE, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list of lists of estimated crossover locations, with
#'     crossovers placed at the midpoint of the intervals that contain
#'     them.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' pr <- calc_genoprob(iron, map, error_prob=0.002, map_function="c-f")
#' g <- maxmarg(pr)
#' pos <- locate_xo(g, iron$gmap)
#'
#' @seealso [count_xo()]
#'
#' @export

locate_xo <-
    function(geno, map, quiet=TRUE, cores=1)
{
    if(is.null(geno)) stop("geno is NULL")
    if(is.null(map)) stop("map is NULL")

    crosstype <- attr(geno, "crosstype")
    if(is.null(crosstype))
        stop("Input geno needs to include a crosstype attribute.")
    is_x_chr <- attr(geno, "is_x_chr")
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(geno))
    names(is_x_chr) <- names(geno)

    if(length(geno) != length(map) ||
       names(geno) != names(map)) { # force matching chromosomes
        chr <- find_common_ids(names(geno), names(map))
        if(length(chr)==0)
            stop("geno and map have no chromosomes in common")
        geno <- geno[,chr]
        map <- map[chr]
    }
    for(i in seq_along(geno)) { # get matching markers
        if(ncol(geno[[i]]) != length(map[[i]]) ||
           !all(colnames(geno[[i]]) == names(map[[i]]))) {
            mar <- find_common_ids(colnames(geno[[i]]), names(map[[i]]))
            if(length(mar) == 0) {
                warning("No markers in common on chr ", names(map)[i])
                geno <- geno[-i]
                map <- map[-i]
            }
            else {
                geno[[i]] <- geno[[i]][,mar,drop=FALSE]
                map[[i]] <- map[[i]][mar]
            }
        }
    }
    if(length(geno) == 0)
        stop("geno and map have no chromosomes/markers in common")
    is_x_chr <- is_x_chr[names(geno)]

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    by_chr_func <- function(chr) {
        result <- .locate_xo(t(geno[[chr]]), map[[chr]], crosstype, is_x_chr[chr])
        names(result) <- rownames(geno[[chr]])
        result
    }

    result <- cluster_lapply(cores, seq(along=geno), by_chr_func)
    names(result) <- names(geno)
    result
}
