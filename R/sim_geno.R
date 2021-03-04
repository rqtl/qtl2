# sim_geno
#' Simulate genotypes given observed marker data
#'
#' Uses a hidden Markov model to simulate from the joint distribution
#' Pr(g | O) where g is the underlying sequence of true genotypes and
#' O is the observed multipoint marker data, with possible allowance
#' for genotyping errors.
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param map Genetic map of markers. May include pseudomarker
#' locations (that is, locations that are not within the marker
#' genotype data). If NULL, the genetic map in `cross` is used.
#' @param n_draws Number of simulations to perform.
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#' @param lowmem If `FALSE`, split individuals into groups with
#' common sex and crossinfo and then precalculate the transition
#' matrices for a chromosome; potentially a lot faster but using more
#' memory.
#' @param quiet If `FALSE`, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return An object of class `"sim_geno"`: a list of three-dimensional arrays of imputed genotypes,
#' individuals x positions x draws. Also contains three attributes:
#' * `crosstype` - The cross type of the input `cross`.
#' * `is_x_chr` - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input `cross`.
#' * `alleles` - Vector of allele codes, from input
#'     `cross`.
#'
#' @details
#'  After performing the backward equations, we draw from
#'  \eqn{Pr(g_1 = v | O)}{Pr(g[1] = v | O)} and then \eqn{Pr(g_{k+1} = v |
#'    O, g_k = u)}{Pr(g[k+1] = v | O, g[k] = u)}.
#'
#' @seealso [cbind.sim_geno()], [rbind.sim_geno()]
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map_w_pmar <- insert_pseudomarkers(grav2$gmap, step=1)
#' draws <- sim_geno(grav2, map_w_pmar, n_draws=4, error_prob=0.002)

sim_geno <-
function(cross, map=NULL, n_draws=1, error_prob=1e-4,
         map_function=c("haldane", "kosambi", "c-f", "morgan"),
         lowmem=FALSE, quiet=TRUE, cores=1)
{
    # check inputs
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)

    if(!is_pos_number(n_draws)) stop("n_draws should be a single positive integer")
    if(!is_nonneg_number(error_prob)) stop("error_prob should be a single non-negative number")

    if(!lowmem)
        return(sim_geno2(cross=cross, map=map, n_draws=n_draws,
                         error_prob=error_prob, map_function=map_function, quiet=quiet,
                         cores=cores))

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # no more messages
    }

    # pseudomarker map
    if(is.null(map)) {
        if(is.null(cross$gmap)) stop("If cross does not contain a genetic map, map must be provided.")
        map <- insert_pseudomarkers(cross$gmap)
    }
    # possibly subset the map
    if(length(map) != length(cross$geno) || !all(names(map) == names(cross$geno))) {
        chr <- names(cross$geno)
        if(!all(chr %in% names(map)))
            stop("map doesn't contain all of the necessary chromosomes")
        map <- map[chr]
    }
    # calculate marker index object
    index <- create_marker_index(lapply(cross$geno, colnames), map)

    rf <- map2rf(map, map_function)

    # deal with missing information
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_group_func <- function(i) {
        dr <- .sim_geno(cross$crosstype, t(cross$geno[[chr]][group[[i]],,drop=FALSE]),
                        founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[group[[i]]],
                        t(cross$cross_info[group[[i]],,drop=FALSE]), rf[[chr]], index[[chr]],
                        error_prob, n_draws)
        aperm(dr, c(3,1,2))
    }

    group <- parallel::splitIndices(nrow(cross$geno[[1]]), n_cores(cores))
    groupindex <- seq(along=group)

    draws <- vector("list", length(cross$geno))
    names(draws) <- names(cross$geno)
    for(chr in seq(along=cross$geno)) {
        if(!quiet) message("Chr ", names(cross$geno)[chr])

        if(n_cores(cores)==1) { # no parallel processing
            # calculations in one group
            draws[[chr]] <- by_group_func(1)
        }
        else {
            # calculations in parallel
            temp <- cluster_lapply(cores, groupindex, by_group_func)

            # paste them back together
            d <- vapply(temp, dim, rep(0,3))
            nr <- sum(d[1,])
            draws[[chr]] <- array(dim=c(nr, d[2,1], d[3,1]))
            for(i in groupindex)
                draws[[chr]][group[[i]],,] <- temp[[i]]
        }

        dimnames(draws[[chr]]) <- list(rownames(cross$geno[[chr]]),
                                       names(map[[chr]]),
                                       NULL)

    }

    names(draws) <- names(cross$geno)
    attr(draws, "crosstype") <- cross$crosstype
    attr(draws, "is_x_chr") <- cross$is_x_chr
    attr(draws, "alleles") <- cross$alleles

    class(draws) <- c("sim_geno", "list")
    draws
}
