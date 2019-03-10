# est_map
#' Estimate genetic maps
#'
#' Uses a hidden Markov model to re-estimate the genetic map for an
#' experimental cross, with possible allowance for genotyping errors.
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#' @param lowmem If `FALSE`, precalculate initial and emission
#' probabilities, and at each iteration calculate the transition
#' matrices for a chromosome; potentially a lot faster but using
#' more memory. Needs to be tailored somewhat to cross type. For
#' example, multi-way RIL may need to reorder the transition
#' matrix according to cross order, and AIL and DO need separate
#' transition matrices for each generation.
#' @param maxit Maximum number of iterations in EM algorithm.
#' @param tol Tolerance for determining convergence
#' @param quiet If `FALSE`, print progress messages.
#' @param save_rf If `TRUE`, save the estimated recombination
#' fractions as an attribute (`"rf"`) of the result.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list of numeric vectors, with the estimated marker
#' locations (in cM). The location of the initial marker on each
#' chromosome is kept the same as in the input `cross`.
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
#' \dontshow{grav2 <- grav2[,"3"]}
#' gmap <- est_map(grav2, error_prob=0.002)

est_map <-
function(cross, error_prob=1e-4,
         map_function=c("haldane", "kosambi", "c-f", "morgan"),
         lowmem=FALSE, maxit=10000, tol=1e-6, quiet=TRUE, save_rf=FALSE,
         cores=1)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    map_function <- match.arg(map_function)
    if(!is_nonneg_number(error_prob) || error_prob > 1) stop("error_prob must be a single number in [0,1]")
    if(!is_nonneg_number(maxit)) stop("maxit must be a single non-negative number")
    if(!is_pos_number(tol)) stop("tol must be a single positive number")

    # deal with missing information
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    map <- vector("list", length(cross$gmap))

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # no more messages
    }

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_chr_func <- function(chr) {
        # the following avoids a warning in R CMD check
        . <- "avoid R CMD check warning"

        if(!quiet) message("Chr ", names(cross$geno)[chr])

        gmap <- cross$gmap[[chr]]

        # omit individuals with < 2 genotypes
        geno <- cross$geno[[chr]]
        ntyped <- rowSums(geno>0)
        keep <- ntyped >= 2
        geno <- t(geno[keep,,drop=FALSE])
        is_female <- cross$is_female[keep]
        cross_info <- t(cross$cross_info[keep,,drop=FALSE])

        rf_start <- map2rf(gmap) # positions to inter-marker rec frac
        if(lowmem)
            rf <- .est_map(cross$crosstype, geno, founder_geno[[chr]],
                           is_x_chr[chr], is_female, cross_info,
                           rf_start, error_prob, maxit, tol, !quiet)
        else {
            # groups of individuals with common sex and cross_info
            sex_crossinfo <- paste(is_female, apply(cross_info, 2, paste, collapse=":"), sep=":")
            unique_cross_group <- unique(sex_crossinfo)
            cross_group <- match(sex_crossinfo, unique_cross_group)-1 # indexes start at 0
            unique_cross_group <- match(seq_along(unique_cross_group)-1, cross_group)-1 # again start at 0

            rf <- .est_map2(cross$crosstype, geno, founder_geno[[chr]],
                            is_x_chr[chr], is_female, cross_info,
                            cross_group, unique_cross_group,
                            rf_start, error_prob, maxit, tol, !quiet)
        }

        loglik <- attr(rf, "loglik")
        map <- cumsum(c(gmap[1], imf(rf, map_function))) # rec frac to positions

        names(map) <- names(gmap)
        attr(map, "loglik") <- loglik
        if(save_rf) {
            attr(rf, "loglik") <- NULL
            attr(map, "rf") <- rf
        }

        map
    }

    chrs <- seq(along=map)
    map <- cluster_lapply(cores, chrs, by_chr_func) # if cores==1, uses lapply

    names(map) <- names(cross$gmap)
    attr(map, "is_x_chr") <- is_x_chr

    if(save_rf) { # put est'd rec fracs as single attribute, as list
        rfs <- lapply(map, function(a) attr(a, "rf"))
        names(rfs) <- names(map)
        attr(map, "rf") <- rfs

        # strip off the individual rf attributes
        for(i in seq_along(map))
            attr(map[[i]], "rf") <- NULL
    }

    map
}
