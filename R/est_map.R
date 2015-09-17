# est_map
#' Estimate genetic maps
#'
#' Uses a hidden Markov model to re-estimate the genetic map for an
#' experimental cross, with possible allowance for genotyping errors.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#' @param maxit Maximum number of iterations in EM algorithm.
#' @param tol Tolerance for determining convergence
#' @param quiet If \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A list of numeric vectors, with the estimated marker
#' locations (in cM). The location of the initial marker on each
#' chromosome is kept the same as in the input \code{cross}.
#'
#' @details
#' The map is estimated assuming no crossover interference,
#' but a map function (by default, Haldane's) is used to derive the genetic distances.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' gmap <- est_map(grav2, error_prob=0.002)

est_map <-
function(cross, error_prob=1e-4,
         map_function=c("haldane", "kosambi", "c-f", "morgan"),
         maxit=10000, tol=1e-6, quiet=TRUE,
         cores=1)
{
    map_function <- match.arg(map_function)
    if(error_prob < 0) stop("error_prob must be >= 0")
    if(maxit < 0) stop("maxit must be >= 0")
    if(tol <= 0) stop("tol must be > 0")

    # deal with missing information
    n.ind <- nrow(cross$geno[[1]])
    chrnames <- names(cross$geno)
    cross_info <- handle_null_crossinfo(cross$cross_info, n.ind)
    is_female <- handle_null_isfemale(cross$is_female, n.ind)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)

    cross_info <- t(cross_info)

    map <- vector("list", length(cross$gmap))

    if("cluster" %in% class(cores) && "SOCKcluster" %in% class(cores)) { # cluster already set
        cluster_ready <- TRUE
        if(!quiet) message(" - Using ", length(cores), " cores.")
        quiet <- TRUE # no more messages
    } else {
        cluster_ready <- FALSE
        if(cores==0) cores <- parallel::detectCores() # if 0, detect cores
        if(cores > 1) {
            if(!quiet) message(" - Using ", cores, " cores.")
            quiet <- TRUE # no more messages
        }
    }

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_chr_func <- function(chr) {
        # the following avoids a warning in R CMD check
        . <- "avoid R CMD check warning"

        if(!quiet) cat(paste0("Chr ", names(cross$geno)[chr], ":\n"))

        gmap <- cross$gmap[[chr]]

        # omit individuals with < 2 genotypes
        geno <- cross$geno[[chr]]
        ntyped <- rowSums(geno>0)
        keep <- (ntyped >= 2)

        rf_start <- mf(diff(gmap), map_function) # positions to inter-marker rec frac
        rf <- .est_map(cross$crosstype, t(cross$geno[[chr]][keep,,drop=FALSE]),
                       founder_geno[[chr]], is_x_chr[chr], is_female[keep],
                       cross_info[,keep,drop=FALSE],
                       rf_start,
                       error_prob, maxit, tol, !quiet)

        loglik <- attr(rf, "loglik")
        map <- cumsum(c(gmap[1], imf(rf, map_function))) # rec frac to positions

        names(map) <- names(gmap)
        attr(map, "loglik") <- loglik

        map
    }

    chrs <- seq(along=map)
    if(!cluster_ready && cores<=1) { # no parallel processing
        map <- lapply(chrs, by_chr_func)
    }
    else if(cluster_ready || Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        if(!cluster_ready) {
            cores <- parallel::makeCluster(cores)
            on.exit(parallel::stopCluster(cores))
        }
        map <- parallel::clusterApply(cores, chrs, by_chr_func)
    }
    else {
        map <- parallel::mclapply(chrs, by_chr_func, mc.cores=cores)
    }

    names(map) <- names(cross$gmap)
    attr(map, "is_x_chr") <- is_x_chr
    map
}
