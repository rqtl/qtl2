# calc_genoprob
#' Calculate conditional genotype probabilities
#'
#' Uses a hidden Markov model to calculate the probabilities of the
#' true underlying genotypes given the observed multipoint marker
#' data, with possible allowance for genotyping errors.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
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
#' @return A list of three-dimensional arrays of probabilities,
#' individuals x genotypes x positions. (Note that the arrangement is
#' different from R/qtl.) The genetic map and cross information are
#' included as attributes.
#'
#' @details
#'   Let \eqn{O_k}{O[k]} denote the observed marker genotype at position
#'  \eqn{k}, and \eqn{g_k}{g[k]} denote the corresponding true underlying
#'  genotype.
#'
#'  We use the forward-backward equations to calculate
#'  \eqn{\alpha_{kv} = \log Pr(O_1, \ldots, O_k, g_k = v)}{%
#'    a[k][v] = log Pr(O[1], \ldots, O[k], g[k] = v)}
#'  and
#'  \eqn{\beta_{kv} = \log Pr(O_{k+1}, \ldots, O_n | g_k = v)}{%
#'    b[k][v] = log Pr(O[k+1], \ldots, O[n] | g[k] = v)}
#'
#'  We then obtain
#'  \eqn{Pr(g_k | O_1, \ldots, O_n) = \exp(\alpha_{kv} + \beta_{kv}) / s}{%
#'    Pr(g[k] | O[1], \ldots, O[n] = exp(a[k][v] + b[k][v]) / s}
#'  where
#'  \eqn{s = \sum_v \exp(\alpha_{kv} + \beta_{kv})}{%
#'    s = sum_v exp(a[k][v] + b[k][v])}
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

calc_genoprob <-
function(cross, step=0, off_end=0, stepwidth=c("fixed", "max"), pseudomarker_map,
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
    rf <- lapply(map, function(m) mf(diff(m), map_function))

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

        pr <- .calc_genoprob(cross$crosstype, t(cross$geno[[chr]]),
                             founder_geno[[chr]], cross$is_x_chr[chr], is_female,
                             cross_info, rf[[chr]], attr(map[[chr]], "index"),
                             error_prob)
        pr <- aperm(pr, c(2,1,3))

        dimnames(pr) <- list(rownames(cross$geno[[chr]]),
                             NULL, # FIX ME: need genotype names in here
                             names(map[[chr]]))
        pr
    }

    chrs <- seq(along=map)
    if(!cluster_ready && cores<=1) { # no parallel processing
        probs <- lapply(chrs, by_chr_func)
    }
    else if(cluster_ready || Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
        if(!cluster_ready) {
            cores <- parallel::makeCluster(cores)
            on.exit(parallel::stopCluster(cores))
        }
        probs <- parallel::clusterApply(cores, chrs, by_chr_func)
    }
    else {
        probs <- parallel::mclapply(chrs, by_chr_func, mc.cores=cores)
    }

    names(probs) <- names(cross$gmap)
    attr(probs, "map") <- map
    attr(probs, "is_x_chr") <- cross$is_x_chr
    attr(probs, "crosstype") <- cross$crosstype
    attr(probs, "sex") <- cross$sex
    attr(probs, "cross_info") <- cross$cross_info

    probs
}

# create empty set of matrices for founder genotype data
create_empty_founder_geno <-
function(geno)
{
    result <- vector("list", length(geno))
    names(result) <- names(geno)
    for(i in seq(along=geno)) {
        result[[i]] <- matrix(0L, nrow=0, ncol=ncol(geno[[i]]))
        colnames(result[[i]]) <- colnames(geno[[i]])
    }
    result
}
