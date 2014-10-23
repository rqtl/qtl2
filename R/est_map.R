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
#' @param quiet If true, don't print any messages (
#' @param n_cores Number of CPU cores to use, for parallel calculations.
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
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' gmap <- est_map(grav2, error_prob=0.002, n_cores=1)

est_map <-
function(cross, error_prob=1e-4,
         map_function=c("haldane", "kosambi", "c-f", "morgan"),
         maxit=10000, tol=1e-6, quiet=TRUE,
         n_cores=1)
{
    map_function <- match.arg(map_function)
    if(error_prob < 0) stop("error_prob must be >= 0")
    if(maxit < 0) stop("maxit must be >= 0")
    if(tol <= 0) stop("tol must be > 0")

    cross_info <- t(cross$cross_info)
    is_female <- cross$is_female

    map <- vector("list", length(cross$gmap))

    if(n_cores > 1) quiet <- TRUE

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

        rf <- .est_map(cross$crosstype, t(cross$geno[[chr]][keep,,drop=FALSE]),
                       founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[keep],
                       cross_info[,keep,drop=FALSE],
                       diff(gmap) %>% mf(map_function), # positions to inter-marker rec frac
                       error_prob, maxit, tol, !quiet)

        loglik <- attr(rf, "loglik")
        map <- imf(rf, map_function) %>% c(gmap[1], .) %>% cumsum() # rec frac to positions

        names(map) <- names(gmap)
        attr(map, "loglik") <- loglik

        map
    }

    chrs <- seq(along=map)
    if(n_cores<=1) { # no parallel processing
        map <- lapply(chrs, by_chr_func)
    }
    else if(Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))
        map <- parallel::clusterApply(cl, chrs, by_chr_func)
    }
    else {
        map <- parallel::mclapply(chrs, by_chr_func, mc.cores=n_cores)
    }

    names(map) <- names(cross$gmap)
    map
}
