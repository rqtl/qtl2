# sim_geno
#' Simulate genotypes given observed marker data
#'
#' Uses a hidden Markov model to simulate from the joint distribution
#' Pr(g | O) where g is the underlying sequence of true genotypes and
#' O is the observed multipoint marker data, with possible allowance
#' for genotyping errors.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param n_draws Number of simulations to perform.
#' @param step Distance between pseudomarkers and markers; if
#' \code{step=0} no pseudomarkers are inserted.
#' @param off_end Distance beyond terminal markers in which to insert
#' pseudomarkers.
#' @param stepwidth Indicates whether to use a fixed grid
#' (\code{stepwidth="fixed"}) or to use the maximal distance between
#' pseudomarkers to ensure that no two adjacent markers/pseudomarkers
#' are more than \code{step} apart.
#' @param pseudomarker_map A map of pseudomarker locations; if provided the
#' \code{step}, \code{off_end}, and \code{stepwidth} arguments are
#' ignored.
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#' @param quiet If \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A list of three-dimensional arrays of imputed genotypes,
#' individuals x positions x draws.  The genetic map and cross
#' information are included as attributes.
#'
#' @details
#'  After performing the backward equations, we draw from
#'  \eqn{Pr(g_1 = v | O)}{Pr(g[1] = v | O)} and then \eqn{Pr(g_{k+1} = v |
#'    O, g_k = u)}{Pr(g[k+1] = v | O, g[k] = u)}.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' draws <- sim_geno(grav2, n_draws=4, step=1, error_prob=0.002)

sim_geno <-
function(cross, n_draws=1, step=0, off_end=0, stepwidth=c("fixed", "max"), pseudomarker_map,
         error_prob=1e-4, map_function=c("haldane", "kosambi", "c-f", "morgan"),
         quiet=TRUE, cores=1)
{
    # check inputs
    if(class(cross) != "cross2")
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)
    stepwidth <- match.arg(stepwidth)

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

    # construct map at which to do the calculations
    if(missing(pseudomarker_map))
        pseudomarker_map <- NULL
    # tolerance for matching marker and pseudomarker positions
    tol <- ifelse(step==0 || step>1, 0.01, step/100)
    # create the combined marker/pseudomarker map
    map <- insert_pseudomarkers(cross$gmap, step, off_end, stepwidth,
                                pseudomarker_map, tol)

    probs <- vector("list", length(map))
    rf <- map2rf(map, map_function)

    # deal with missing information
    n.ind <- nrow(cross$geno[[1]])
    chrnames <- names(cross$geno)
    cross_info <- handle_null_crossinfo(cross$cross_info, n.ind)
    is_female <- handle_null_isfemale(cross$is_female, n.ind)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)

    cross_info <- t(cross$cross_info)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_chr_func <- function(chr) {
        if(!quiet) cat("Chr ", names(cross$geno)[chr], "\n")

        dr <- .sim_geno(cross$crosstype, t(cross$geno[[chr]]),
                        founder_geno[[chr]], cross$is_x_chr[chr], is_female,
                        cross_info, rf[[chr]], attr(map[[chr]], "index"),
                        error_prob, n_draws)
        dr <- aperm(dr, c(3,1,2))

        dimnames(dr) <- list(rownames(cross$geno[[chr]]),
                             names(map[[chr]]),
                             NULL)
        dr
    }

    chrs <- seq(along=map)
    if(!cluster_ready && cores<=1) { # no parallel processing
        draws <- lapply(chrs, by_chr_func)
    }
    else if(cluster_ready || Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        if(!cluster_ready) {
            cores <- parallel::makeCluster(cores)
            on.exit(parallel::stopCluster(cores))
        }
        draws <- parallel::clusterApply(cores, chrs, by_chr_func)
    }
    else {
        draws <- parallel::mclapply(chrs, by_chr_func, mc.cores=cores)
    }

    names(draws) <- names(cross$gmap)
    attr(draws, "map") <- map
    attr(draws, "is_x_chr") <- cross$is_x_chr
    attr(draws, "crosstype") <- cross$crosstype
    attr(draws, "sex") <- cross$sex
    attr(draws, "cross_info") <- cross$cross_info

    class(draws) <- c("sim_geno", "list")
    draws
}
