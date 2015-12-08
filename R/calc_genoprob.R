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

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
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
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_group_func <- function(i) {
        pr <- .calc_genoprob(cross$crosstype, t(cross$geno[[chr]][group[[i]],,drop=FALSE]),
                             founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[group[[i]]],
                             t(cross$cross_info[group[[i]],,drop=FALSE]), rf[[chr]], attr(map[[chr]], "index"),
                             error_prob)
        aperm(pr, c(2,1,3))
    }

    group <- vec4parallel(nrow(cross$geno[[1]]), n_cores(cores))
    groupindex <- seq(along=group)

    probs <- vector("list", length(cross$geno))
    for(chr in seq(along=cross$geno)) {
        if(!quiet) cat("Chr ", names(cross$geno)[chr], "\n")

        if(n_cores(cores) == 1) { # no parallel processing
            # calculations in one group
            probs[[chr]] <- by_group_func(1)
        }
        else {
            # calculations in parallel
            temp <- run_by_cluster(cores, groupindex, by_group_func)

            # paste them back together
            d <- vapply(temp, dim, rep(0,3))
            nr <- sum(d[1,])
            probs[[chr]] <- array(dim=c(nr, d[2,1], d[3,1]))
            for(i in groupindex)
                probs[[chr]][group[[i]],,] <- temp[[i]]
        }

        dimnames(probs[[chr]]) <- list(rownames(cross$geno[[chr]]),
                                       NULL, # FIX ME: need genotype names in here
                                       names(map[[chr]]))

    }

    names(probs) <- names(cross$gmap)
    attr(probs, "map") <- map
    attr(probs, "is_x_chr") <- cross$is_x_chr
    attr(probs, "crosstype") <- cross$crosstype
    attr(probs, "sex") <- cross$sex
    attr(probs, "cross_info") <- cross$cross_info

    class(probs) <- c("calc_genoprob", "list")

    probs
}
